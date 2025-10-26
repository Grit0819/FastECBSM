package btcec

// References:
//   [SECG]: Recommended Elliptic Curve Domain Parameters
//     http://www.secg.org/sec2-v2.pdf
//
//   [GECC]: Guide to Elliptic Curve Cryptography (Hankerson, Menezes, Vanstone)

// This package operates, internally, on Jacobian coordinates. For a given
// (x, y) position on the curve, the Jacobian coordinates are (x1, y1, z1)
// where x = x1/z1² and y = y1/z1³. The greatest speedups come when the whole
// calculation can be performed within the transform (as in ScalarMult and
// ScalarBaseMult). But even for Add and Double, it's faster to apply and
// reverse the transform than to operate in affine coordinates.

import (
	"crypto/elliptic"
	"fmt"
	"math/big"
	"sync"
	"time"
)

var (
	// fieldOne is simply the integer 1 in field representation.  It is
	// used to avoid needing to create it multiple times during the internal
	// arithmetic.
	fieldOne  = new(FieldVal).SetInt(1)
	fieldZero = new(FieldVal).SetInt(0)
)

// KoblitzCurve supports a koblitz curve implementation that fits the ECC Curve
// interface from crypto/elliptic.
type KoblitzCurve struct {
	*elliptic.CurveParams

	// q is the value (P+1)/4 used to compute the square root of field
	// elements.
	q *big.Int

	H         int      // cofactor of the curve.
	halfOrder *big.Int // half the order N

	// fieldB is the constant B of the curve as a FieldVal.
	fieldB *FieldVal

	// byteSize is simply the bit size / 8 and is provided for convenience
	// since it is calculated repeatedly.
	byteSize int

	// bytePoints
	bytePoints *[32][256][3]FieldVal

	// The next 6 values are used specifically for endomorphism
	// optimizations in ScalarMult.

	// lambda must fulfill lambda^3 = 1 mod N where N is the order of G.
	lambda *big.Int

	// beta must fulfill beta^3 = 1 mod P where P is the prime field of the
	// curve.
	beta *FieldVal

	// See the EndomorphismVectors in gensecp256k1.go to see how these are
	// derived.
	a1 *big.Int
	b1 *big.Int
	a2 *big.Int
	b2 *big.Int
}

// Params returns the parameters for the curve.
func (curve *KoblitzCurve) Params() *elliptic.CurveParams {
	return curve.CurveParams
}

// bigAffineToField takes an affine point (x, y) as big integers and converts
// it to an affine point as field values.

func (curve *KoblitzCurve) bigAffineToField(x, y *big.Int) (*FieldVal, *FieldVal) {
	x3, y3 := new(FieldVal), new(FieldVal)
	x3.SetByteSlice(x.Bytes())
	y3.SetByteSlice(y.Bytes())
	return x3, y3
}

func (curve *KoblitzCurve) FbigAffineToField(x, y *big.Int) (*FieldVal, *FieldVal) {
	x3, y3 := new(FieldVal), new(FieldVal)
	x3.SetByteSlice(x.Bytes())
	y3.SetByteSlice(y.Bytes())
	return x3, y3
}

func (curve *KoblitzCurve) bigAffineToField1(x, y *big.Int) (*FieldVal, *FieldVal, *FieldVal) {
	x3, y3, z3 := new(FieldVal), new(FieldVal), new(FieldVal)
	x3.SetByteSlice(x.Bytes())
	y3.SetByteSlice(y.Bytes())
	z3.SetInt(1)
	return x3, y3, z3
}

func (curve *KoblitzCurve) BigAffineToField(x, y *big.Int) (*FieldVal, *FieldVal, *FieldVal) {
	x3, y3, z3 := new(FieldVal), new(FieldVal), new(FieldVal)
	x3.SetByteSlice(x.Bytes())
	y3.SetByteSlice(y.Bytes())
	z3.SetInt(1)
	return x3, y3, z3
}

// 自己写的
func (curve *KoblitzCurve) FieldToBigAffine(x, y *FieldVal) (*big.Int, *big.Int) {
	// Inversions are expensive and both point addition and point doubling
	// are faster when working with points that have a z value of one.  So,
	// if the point needs to be converted to affine, go ahead and normalize
	// the point itself at the same time as the calculation is the same.

	// Normalize the x and y values.
	x.Normalize()
	y.Normalize()

	// Convert the field values for the now affine point to big.Ints.
	x3, y3 := new(big.Int), new(big.Int)
	x3.SetBytes(x.Bytes()[:])
	y3.SetBytes(y.Bytes()[:])
	return x3, y3
}

// fieldJacobianToBigAffine takes a Jacobian point (x, y, z) as field values and
// converts it to an affine point as big integers.

func (curve *KoblitzCurve) fieldJacobianToBigAffine(x, y, z *FieldVal) (*big.Int, *big.Int) {
	// Inversions are expensive and both point addition and point doubling
	// are faster when working with points that have a z value of one.  So,
	// if the point needs to be converted to affine, go ahead and normalize
	// the point itself at the same time as the calculation is the same.
	var zInv, tempZ FieldVal
	zInv.Set(z).Inverse()   // zInv = Z^-1
	tempZ.SquareVal(&zInv)  // tempZ = Z^-2
	x.Mul(&tempZ)           // X = X/Z^2 (mag: 1)
	y.Mul(tempZ.Mul(&zInv)) // Y = Y/Z^3 (mag: 1)
	z.SetInt(1)             // Z = 1 (mag: 1)

	// Normalize the x and y values.
	x.Normalize()
	y.Normalize()

	// Convert the field values for the now affine point to big.Ints.
	x3, y3 := new(big.Int), new(big.Int)
	x3.SetBytes(x.Bytes()[:])
	y3.SetBytes(y.Bytes()[:])
	return x3, y3
}

// 利用逆树来进行还原
func (curve *KoblitzCurve) fieldJacobianToBigAffineBatch(x, y []*FieldVal, ZinvTree [Treelen]*FieldVal) ([]*big.Int, []*big.Int) {
	// Inversions are expensive and both point addition and point doubling
	// are faster when working with points that have a z value of one.  So,
	// if the point needs to be converted to affine, go ahead and normalize
	// the point itself at the same time as the calculation is the same.
	var resx, resy []*big.Int
	//批量进行还原
	for i := 0; i < len(x); i++ {
		var tempZ FieldVal
		// zInv = Z^-1
		tempZ.SquareVal(ZinvTree[i])     // tempZ = Z^-2
		x[i].Mul(&tempZ)                 // X = X/Z^2 (mag: 1)
		y[i].Mul(tempZ.Mul(ZinvTree[i])) // Y = Y/Z^3 (mag: 1)

		// Normalize the x and y values.
		x[i].Normalize()
		y[i].Normalize()

		// Convert the field values for the now affine point to big.Ints.
		resx = append(resx, new(big.Int))
		resy = append(resy, new(big.Int))
		resx[i].SetBytes(x[i].Bytes()[:])
		resy[i].SetBytes(y[i].Bytes()[:])

	}

	return resx, resy

}

func (curve *KoblitzCurve) FieldJacobianToBigAffine(x, y, z *FieldVal) (*big.Int, *big.Int) {
	// Inversions are expensive and both point addition and point doubling
	// are faster when working with points that have a z value of one.  So,
	// if the point needs to be converted to affine, go ahead and normalize
	// the point itself at the same time as the calculation is the same.
	var zInv, tempZ FieldVal
	zInv.Set(z).Inverse()   // zInv = Z^-1
	tempZ.SquareVal(&zInv)  // tempZ = Z^-2
	x.Mul(&tempZ)           // X = X/Z^2 (mag: 1)
	y.Mul(tempZ.Mul(&zInv)) // Y = Y/Z^3 (mag: 1)
	z.SetInt(1)             // Z = 1 (mag: 1)

	// Normalize the x and y values.
	x.Normalize()
	y.Normalize()

	// Convert the field values for the now affine point to big.Ints.
	x3, y3 := new(big.Int), new(big.Int)
	x3.SetBytes(x.Bytes()[:])
	y3.SetBytes(y.Bytes()[:])
	return x3, y3
}

func (curve *KoblitzCurve) fieldJacobianToBigAffineWithoutInv(x, y, inv_z *FieldVal) (*big.Int, *big.Int) {
	// Inversions are expensive and both point addition and point doubling
	// are faster when working with points that have a z value of one.  So,
	// if the point needs to be converted to affine, go ahead and normalize
	// the point itself at the same time as the calculation is the same.
	var zInv, tempZ FieldVal
	zInv.Set(inv_z)
	tempZ.SquareVal(&zInv)  // tempZ = Z^-2
	x.Mul(&tempZ)           // X = X/Z^2 (mag: 1)
	y.Mul(tempZ.Mul(&zInv)) // Y = Y/Z^3 (mag: 1)
	// Normalize the x and y values.
	x.Normalize()
	y.Normalize()

	// Convert the field values for the now affine point to big.Ints.
	x3, y3 := new(big.Int), new(big.Int)
	x3.SetBytes(x.Bytes()[:])
	y3.SetBytes(y.Bytes()[:])
	return x3, y3
}

func (curve *KoblitzCurve) fieldJacobianToBigAffineXWithoutInv(x, inv_z *FieldVal) *big.Int {
	// Inversions are expensive and both point addition and point doubling
	// are faster when working with points that have a z value of one.  So,
	// if the point needs to be converted to affine, go ahead and normalize
	// the point itself at the same time as the calculation is the same.
	var zInv, tempZ FieldVal
	zInv.Set(inv_z)
	tempZ.SquareVal(&zInv) // tempZ = Z^-2
	x.Mul(&tempZ)          // X = X/Z^2 (mag: 1)
	x.Normalize()
	// Convert the field values for the now affine point to big.Ints.
	x3 := new(big.Int)
	x3.SetBytes(x.Bytes()[:])
	return x3
}

// IsOnCurve returns boolean if the point (x,y) is on the curve.
// Part of the elliptic.Curve interface. This function differs from the
// crypto/elliptic algorithm since a = 0 not -3.
func (curve *KoblitzCurve) IsOnCurve(x, y *big.Int) bool {
	// Convert big ints to field values for faster arithmetic.
	fx, fy := curve.bigAffineToField(x, y)

	// Elliptic curve equation for secp256k1 is: y^2 = x^3 + 7
	y2 := new(FieldVal).SquareVal(fy).Normalize()
	result := new(FieldVal).SquareVal(fx).Mul(fx).AddInt(7).Normalize()
	return y2.Equals(result)
}

// Getz3是一个求(x1,y1,1)+(x2,y2,1)=(x3,y3,z3)时只求z3的函数
func (curve *KoblitzCurve) Getz3(x1, x2 *FieldVal) *FieldVal {
	var z3 FieldVal
	z3.Set(x1).Negate(1).Add(x2) // H = X2-X1 (mag: 3)// 设置为x1的取反 然后加x2 == x2 -x1 H
	return &z3                   //直接返回H
	//z3.Set(&h).MulInt(2) // Z3 = 2*H (mag: 6)
	//z3.Normalize()
}

func (curve *KoblitzCurve) Getx3y3(x1, y1, x2, y2, x3, y3 *FieldVal) {
	//x1.Normalize()
	//y1.Normalize()
	//x2.Normalize()
	//y2.Normalize()
	// Calculate X3, Y3, and Z3 according to the intermediate elements
	// breakdown above.
	var h, i, j, r, v FieldVal
	var negJ, neg2V, negX3 FieldVal
	h.Set(x1).Negate(1).Add(x2)                // H = X2-X1 (mag: 3)
	i.SquareVal(&h).MulInt(4)                  // I = 4*H^2 (mag: 4)
	j.Mul2(&h, &i)                             // J = H*I (mag: 1)
	r.Set(y1).Negate(1).Add(y2).MulInt(2)      // r = 2*(Y2-Y1) (mag: 6)
	v.Mul2(x1, &i)                             // V = X1*I (mag: 1)
	negJ.Set(&j).Negate(1)                     // negJ = -J (mag: 2)
	neg2V.Set(&v).MulInt(2).Negate(2)          // neg2V = -(2*V) (mag: 3)
	x3.Set(&r).Square().Add(&negJ).Add(&neg2V) // X3 = r^2-J-2*V (mag: 6)
	negX3.Set(x3).Negate(6)                    // negX3 = -X3 (mag: 7)
	j.Mul(y1).MulInt(2).Negate(2)              // J = -(2*Y1*J) (mag: 3)
	y3.Set(&v).Add(&negX3).Mul(&r).Add(&j)     // Y3 = r*(V-X3)-2*Y1*J (mag: 4)
	// Normalize the resulting field values to a magnitude of 1 as needed.
	x3.Normalize()
	y3.Normalize()
}

func (curve *KoblitzCurve) Getx3(x1, y1, x2, y2, zinv *FieldVal) *big.Int {
	//x1.Normalize()
	//y1.Normalize()
	//x2.Normalize()
	//y2.Normalize()
	var h, i, r, v, x3 FieldVal
	r.Set(y1).Negate(1).Add(y2) // r = (y2-y1)
	h.Mul2(&r, zinv)            // h = (y2-y1)*(x2-x1)^-1
	i.SquareVal(&h)             // I = H^2 (mag: 4)
	v.Set(x1).Add(x2).Negate(1) // V = -(X1+X2) (mag: 1)
	x3.Set(&i).Add(&v)          // X3 = r^2-v (mag: 6)
	x3.Normalize()
	affx3 := new(big.Int)
	affx3.SetBytes(x3.Bytes()[:])
	return affx3
}

/*
func (curve *KoblitzCurve) NewGetx3(x1, y1, x2, y2, zinv, p *FieldVal) (*big.Int, *big.Int) {

	var h, i, r, v, x3, invy2, r1, invy1, h1, i1, invx3 FieldVal
	invy2.Set(y2).Negate(1).Add(p)
	invy1.Set(y1).Negate(1)
	r.Set(&invy1).Add(y2) // r = (y2-y1)
	r1.Set(&invy1).Add(&invy2)
	h.Mul2(&r, zinv) // h = (y2-y1)*(x2-x1)^-1
	h1.Mul2(&r1, zinv)
	i.SquareVal(&h) // I = H^2 (mag: 4)
	i1.SquareVal(&h1)
	v.Set(x1).Add(x2).Negate(1) // V = -(X1+X2) (mag: 1)
	x3.Set(&i).Add(&v)          // X3 = r^2-v (mag: 6)
	invx3.Set(&i1).Add(&v)
	x3.Normalize()
	invx3.Normalize()
	affx3 := new(big.Int).SetBytes(x3.Bytes()[:])
	invaffx3 := new(big.Int).SetBytes(invx3.Bytes()[:])
	return affx3, invaffx3
}*/

// 实际的利用Montgomery Trick 计算两点 x3 = k^2 - x1 - x2   y3 = k(x1 -x3) - y1   a==b  k = (3x1^2 + a)/2y1   a/= b k = (y2 - y1)/(x2 - x1)
func (curve *KoblitzCurve) NewGetx3(x1, y1, x2, y2, zinv, p *FieldVal) (*big.Int, *big.Int) {
	var h, i, r, v, x3, r1, invy1, h1, i1, invx3 FieldVal
	invy1.Set(y1).Negate(1) //-y1
	r.Set(&invy1).Add(y2)   // r = (y2-y1)
	r1.Set(y1).Add(y2)      //y2-y1   这应该是y1 + y2
	h.Mul2(&r, zinv)        // h = (y2-y1)*(x2-x1)^-1   zinv 是逆树 前面存的叶子结点是x2-x1的逆
	h1.Mul2(&r1, zinv)
	i.SquareVal(&h) // I = H^2 (mag: 4) k的平方
	i1.SquareVal(&h1)
	v.Set(x1).Add(x2).Negate(1) // V = -(X1+X2) (mag: 1) 就是-x1 - x2
	x3.Set(&i).Add(&v)          // X3 = r^2+v (mag: 6) x3
	invx3.Set(&i1).Add(&v)
	x3.Normalize()
	invx3.Normalize()
	affx3 := new(big.Int).SetBytes(x3.Bytes()[:])
	invaffx3 := new(big.Int).SetBytes(invx3.Bytes()[:])
	return affx3, invaffx3
}

func (curve *KoblitzCurve) addZ1AndZ2EqualsOne(x1, y1, z1, x2, y2, x3, y3, z3 *FieldVal) {
	x1.Normalize()
	y1.Normalize()
	x2.Normalize()
	y2.Normalize()
	if x1.Equals(x2) {
		if y1.Equals(y2) {
			// Since x1 == x2 and y1 == y2, point doubling must be
			// done, otherwise the addition would end up dividing
			// by zero.
			curve.doubleJacobian(x1, y1, z1, x3, y3, z3)
			return
		}

		// Since x1 == x2 and y1 == -y2, the sum is the point at
		// infinity per the group law.
		x3.SetInt(0)
		y3.SetInt(0)
		z3.SetInt(0)
		return
	}

	// Calculate X3, Y3, and Z3 according to the intermediate elements
	// breakdown above.
	var h, i, j, r, v FieldVal
	var negJ, neg2V, negX3 FieldVal
	h.Set(x1).Negate(1).Add(x2)                // H = X2-X1 (mag: 3)
	i.SquareVal(&h).MulInt(4)                  // I = 4*H^2 (mag: 4)
	j.Mul2(&h, &i)                             // J = H*I (mag: 1)
	r.Set(y1).Negate(1).Add(y2).MulInt(2)      // r = 2*(Y2-Y1) (mag: 6)
	v.Mul2(x1, &i)                             // V = X1*I (mag: 1)
	negJ.Set(&j).Negate(1)                     // negJ = -J (mag: 2)
	neg2V.Set(&v).MulInt(2).Negate(2)          // neg2V = -(2*V) (mag: 3)
	x3.Set(&r).Square().Add(&negJ).Add(&neg2V) // X3 = r^2-J-2*V (mag: 6)
	negX3.Set(x3).Negate(6)                    // negX3 = -X3 (mag: 7)
	j.Mul(y1).MulInt(2).Negate(2)              // J = -(2*Y1*J) (mag: 3)
	y3.Set(&v).Add(&negX3).Mul(&r).Add(&j)     // Y3 = r*(V-X3)-2*Y1*J (mag: 4)
	z3.Set(&h).MulInt(2)                       // Z3 = 2*H (mag: 6)

	// Normalize the resulting field values to a magnitude of 1 as needed.
	x3.Normalize()
	y3.Normalize()
	z3.Normalize()
}

// addZ1EqualsZ2 adds two Jacobian points that are already known to have the
// same z value and stores the result in (x3, y3, z3).  That is to say
// (x1, y1, z1) + (x2, y2, z1) = (x3, y3, z3).  It performs faster addition than
// the generic add routine since less arithmetic is needed due to the known
// equivalence.
func (curve *KoblitzCurve) addZ1EqualsZ2(x1, y1, z1, x2, y2, x3, y3, z3 *FieldVal) {
	// To compute the point addition efficiently, this implementation splits
	// the equation into intermediate elements which are used to minimize
	// the number of field multiplications using a slightly modified version
	// of the method shown at:
	// http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-mmadd-2007-bl
	//
	// In particular it performs the calculations using the following:
	// A = X2-X1, B = A^2, C=Y2-Y1, D = C^2, E = X1*B, F = X2*B
	// X3 = D-E-F, Y3 = C*(E-X3)-Y1*(F-E), Z3 = Z1*A
	//
	// This results in a cost of 5 field multiplications, 2 field squarings,
	// 9 field additions, and 0 integer multiplications.

	// When the x coordinates are the same for two points on the curve, the
	// y coordinates either must be the same, in which case it is point
	// doubling, or they are opposite and the result is the point at
	// infinity per the group law for elliptic curve cryptography.
	x1.Normalize()
	y1.Normalize()
	x2.Normalize()
	y2.Normalize()
	if x1.Equals(x2) {
		if y1.Equals(y2) {
			// Since x1 == x2 and y1 == y2, point doubling must be
			// done, otherwise the addition would end up dividing
			// by zero.
			curve.doubleJacobian(x1, y1, z1, x3, y3, z3)
			return
		}

		// Since x1 == x2 and y1 == -y2, the sum is the point at
		// infinity per the group law.
		x3.SetInt(0)
		y3.SetInt(0)
		z3.SetInt(0)
		return
	}

	// Calculate X3, Y3, and Z3 according to the intermediate elements
	// breakdown above.
	var a, b, c, d, e, f FieldVal
	var negX1, negY1, negE, negX3 FieldVal
	negX1.Set(x1).Negate(1)                // negX1 = -X1 (mag: 2)
	negY1.Set(y1).Negate(1)                // negY1 = -Y1 (mag: 2)
	a.Set(&negX1).Add(x2)                  // A = X2-X1 (mag: 3)
	b.SquareVal(&a)                        // B = A^2 (mag: 1)
	c.Set(&negY1).Add(y2)                  // C = Y2-Y1 (mag: 3)
	d.SquareVal(&c)                        // D = C^2 (mag: 1)
	e.Mul2(x1, &b)                         // E = X1*B (mag: 1)
	negE.Set(&e).Negate(1)                 // negE = -E (mag: 2)
	f.Mul2(x2, &b)                         // F = X2*B (mag: 1)
	x3.Add2(&e, &f).Negate(3).Add(&d)      // X3 = D-E-F (mag: 5)
	negX3.Set(x3).Negate(5).Normalize()    // negX3 = -X3 (mag: 1)
	y3.Set(y1).Mul(f.Add(&negE)).Negate(3) // Y3 = -(Y1*(F-E)) (mag: 4)
	y3.Add(e.Add(&negX3).Mul(&c))          // Y3 = C*(E-X3)+Y3 (mag: 5)
	z3.Mul2(z1, &a)                        // Z3 = Z1*A (mag: 1)

	// Normalize the resulting field values to a magnitude of 1 as needed.
	x3.Normalize()
	y3.Normalize()
}

// addZ2EqualsOne adds two Jacobian points when the second point is already
// known to have a z value of 1 (and the z value for the first point is not 1)
// and stores the result in (x3, y3, z3).  That is to say (x1, y1, z1) +
// (x2, y2, 1) = (x3, y3, z3).  It performs faster addition than the generic
// add routine since less arithmetic is needed due to the ability to avoid
// multiplications by the second point's z value.
func (curve *KoblitzCurve) addZ2EqualsOne(x1, y1, z1, x2, y2, x3, y3, z3 *FieldVal) {
	// To compute the point addition efficiently, this implementation splits
	// the equation into intermediate elements which are used to minimize
	// the number of field multiplications using the method shown at:
	// http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-madd-2007-bl
	//
	// In particular it performs the calculations using the following:
	// Z1Z1 = Z1^2, U2 = X2*Z1Z1, S2 = Y2*Z1*Z1Z1, H = U2-X1, HH = H^2,
	// I = 4*HH, J = H*I, r = 2*(S2-Y1), V = X1*I
	// X3 = r^2-J-2*V, Y3 = r*(V-X3)-2*Y1*J, Z3 = (Z1+H)^2-Z1Z1-HH
	//
	// This results in a cost of 7 field multiplications, 4 field squarings,
	// 9 field additions, and 4 integer multiplications.

	// When the x coordinates are the same for two points on the curve, the
	// y coordinates either must be the same, in which case it is point
	// doubling, or they are opposite and the result is the point at
	// infinity per the group law for elliptic curve cryptography.  Since
	// any number of Jacobian coordinates can represent the same affine
	// point, the x and y values need to be converted to like terms.  Due to
	// the assumption made for this function that the second point has a z
	// value of 1 (z2=1), the first point is already "converted".
	var z1z1, u2, s2 FieldVal
	x1.Normalize()
	y1.Normalize()
	z1z1.SquareVal(z1)                        // Z1Z1 = Z1^2 (mag: 1)
	u2.Set(x2).Mul(&z1z1).Normalize()         // U2 = X2*Z1Z1 (mag: 1)
	s2.Set(y2).Mul(&z1z1).Mul(z1).Normalize() // S2 = Y2*Z1*Z1Z1 (mag: 1)
	if x1.Equals(&u2) {
		if y1.Equals(&s2) {
			// Since x1 == x2 and y1 == y2, point doubling must be
			// done, otherwise the addition would end up dividing
			// by zero.
			curve.doubleJacobian(x1, y1, z1, x3, y3, z3)
			return
		}

		// Since x1 == x2 and y1 == -y2, the sum is the point at
		// infinity per the group law.
		x3.SetInt(0)
		y3.SetInt(0)
		z3.SetInt(0)
		return
	}

	// Calculate X3, Y3, and Z3 according to the intermediate elements
	// breakdown above.
	var h, hh, i, j, r, rr, v FieldVal
	var negX1, negY1, negX3 FieldVal
	negX1.Set(x1).Negate(1)                // negX1 = -X1 (mag: 2)
	h.Add2(&u2, &negX1)                    // H = U2-X1 (mag: 3)
	hh.SquareVal(&h)                       // HH = H^2 (mag: 1)
	i.Set(&hh).MulInt(4)                   // I = 4 * HH (mag: 4)
	j.Mul2(&h, &i)                         // J = H*I (mag: 1)
	negY1.Set(y1).Negate(1)                // negY1 = -Y1 (mag: 2)
	r.Set(&s2).Add(&negY1).MulInt(2)       // r = 2*(S2-Y1) (mag: 6)
	rr.SquareVal(&r)                       // rr = r^2 (mag: 1)
	v.Mul2(x1, &i)                         // V = X1*I (mag: 1)
	x3.Set(&v).MulInt(2).Add(&j).Negate(3) // X3 = -(J+2*V) (mag: 4)
	x3.Add(&rr)                            // X3 = r^2+X3 (mag: 5)
	negX3.Set(x3).Negate(5)                // negX3 = -X3 (mag: 6)
	y3.Set(y1).Mul(&j).MulInt(2).Negate(2) // Y3 = -(2*Y1*J) (mag: 3)
	y3.Add(v.Add(&negX3).Mul(&r))          // Y3 = r*(V-X3)+Y3 (mag: 4)
	z3.Add2(z1, &h).Square()               // Z3 = (Z1+H)^2 (mag: 1)
	z3.Add(z1z1.Add(&hh).Negate(2))        // Z3 = Z3-(Z1Z1+HH) (mag: 4)

	// Normalize the resulting field values to a magnitude of 1 as needed.
	x3.Normalize()
	y3.Normalize()
	z3.Normalize()
}

// addGeneric adds two Jacobian points (x1, y1, z1) and (x2, y2, z2) without any
// assumptions about the z values of the two points and stores the result in
// (x3, y3, z3).  That is to say (x1, y1, z1) + (x2, y2, z2) = (x3, y3, z3).  It
// is the slowest of the add routines due to requiring the most arithmetic.
func (curve *KoblitzCurve) addGeneric(x1, y1, z1, x2, y2, z2, x3, y3, z3 *FieldVal) {
	// To compute the point addition efficiently, this implementation splits
	// the equation into intermediate elements which are used to minimize
	// the number of field multiplications using the method shown at:
	// http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl
	//
	// In particular it performs the calculations using the following:
	// Z1Z1 = Z1^2, Z2Z2 = Z2^2, U1 = X1*Z2Z2, U2 = X2*Z1Z1, S1 = Y1*Z2*Z2Z2
	// S2 = Y2*Z1*Z1Z1, H = U2-U1, I = (2*H)^2, J = H*I, r = 2*(S2-S1)
	// V = U1*I
	// X3 = r^2-J-2*V, Y3 = r*(V-X3)-2*S1*J, Z3 = ((Z1+Z2)^2-Z1Z1-Z2Z2)*H
	//
	// This results in a cost of 11 field multiplications, 5 field squarings,
	// 9 field additions, and 4 integer multiplications.

	// When the x coordinates are the same for two points on the curve, the
	// y coordinates either must be the same, in which case it is point
	// doubling, or they are opposite and the result is the point at
	// infinity.  Since any number of Jacobian coordinates can represent the
	// same affine point, the x and y values need to be converted to like
	// terms.
	var z1z1, z2z2, u1, u2, s1, s2 FieldVal
	z1z1.SquareVal(z1)                        // Z1Z1 = Z1^2 (mag: 1)
	z2z2.SquareVal(z2)                        // Z2Z2 = Z2^2 (mag: 1)
	u1.Set(x1).Mul(&z2z2).Normalize()         // U1 = X1*Z2Z2 (mag: 1)
	u2.Set(x2).Mul(&z1z1).Normalize()         // U2 = X2*Z1Z1 (mag: 1)
	s1.Set(y1).Mul(&z2z2).Mul(z2).Normalize() // S1 = Y1*Z2*Z2Z2 (mag: 1)
	s2.Set(y2).Mul(&z1z1).Mul(z1).Normalize() // S2 = Y2*Z1*Z1Z1 (mag: 1)
	if u1.Equals(&u2) {
		if s1.Equals(&s2) {
			// Since x1 == x2 and y1 == y2, point doubling must be
			// done, otherwise the addition would end up dividing
			// by zero.
			curve.doubleJacobian(x1, y1, z1, x3, y3, z3)
			return
		}

		// Since x1 == x2 and y1 == -y2, the sum is the point at
		// infinity per the group law.
		x3.SetInt(0)
		y3.SetInt(0)
		z3.SetInt(0)
		return
	}

	// Calculate X3, Y3, and Z3 according to the intermediate elements
	// breakdown above.
	var h, i, j, r, rr, v FieldVal
	var negU1, negS1, negX3 FieldVal
	negU1.Set(&u1).Negate(1)               // negU1 = -U1 (mag: 2)
	h.Add2(&u2, &negU1)                    // H = U2-U1 (mag: 3)
	i.Set(&h).MulInt(2).Square()           // I = (2*H)^2 (mag: 2)
	j.Mul2(&h, &i)                         // J = H*I (mag: 1)
	negS1.Set(&s1).Negate(1)               // negS1 = -S1 (mag: 2)
	r.Set(&s2).Add(&negS1).MulInt(2)       // r = 2*(S2-S1) (mag: 6)
	rr.SquareVal(&r)                       // rr = r^2 (mag: 1)
	v.Mul2(&u1, &i)                        // V = U1*I (mag: 1)
	x3.Set(&v).MulInt(2).Add(&j).Negate(3) // X3 = -(J+2*V) (mag: 4)
	x3.Add(&rr)                            // X3 = r^2+X3 (mag: 5)
	negX3.Set(x3).Negate(5)                // negX3 = -X3 (mag: 6)
	y3.Mul2(&s1, &j).MulInt(2).Negate(2)   // Y3 = -(2*S1*J) (mag: 3)
	y3.Add(v.Add(&negX3).Mul(&r))          // Y3 = r*(V-X3)+Y3 (mag: 4)
	z3.Add2(z1, z2).Square()               // Z3 = (Z1+Z2)^2 (mag: 1)
	z3.Add(z1z1.Add(&z2z2).Negate(2))      // Z3 = Z3-(Z1Z1+Z2Z2) (mag: 4)
	z3.Mul(&h)                             // Z3 = Z3*H (mag: 1)

	// Normalize the resulting field values to a magnitude of 1 as needed.
	x3.Normalize()
	y3.Normalize()
}

func (curve *KoblitzCurve) AddGeneric(x1, y1, z1, x2, y2, z2, x3, y3, z3 *FieldVal) {
	// To compute the point addition efficiently, this implementation splits
	// the equation into intermediate elements which are used to minimize
	// the number of field multiplications using the method shown at:
	// http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl
	//
	// In particular it performs the calculations using the following:
	// Z1Z1 = Z1^2, Z2Z2 = Z2^2, U1 = X1*Z2Z2, U2 = X2*Z1Z1, S1 = Y1*Z2*Z2Z2
	// S2 = Y2*Z1*Z1Z1, H = U2-U1, I = (2*H)^2, J = H*I, r = 2*(S2-S1)
	// V = U1*I
	// X3 = r^2-J-2*V, Y3 = r*(V-X3)-2*S1*J, Z3 = ((Z1+Z2)^2-Z1Z1-Z2Z2)*H
	//
	// This results in a cost of 11 field multiplications, 5 field squarings,
	// 9 field additions, and 4 integer multiplications.

	// When the x coordinates are the same for two points on the curve, the
	// y coordinates either must be the same, in which case it is point
	// doubling, or they are opposite and the result is the point at
	// infinity.  Since any number of Jacobian coordinates can represent the
	// same affine point, the x and y values need to be converted to like
	// terms.
	var z1z1, z2z2, u1, u2, s1, s2 FieldVal
	z1z1.SquareVal(z1)                        // Z1Z1 = Z1^2 (mag: 1)
	z2z2.SquareVal(z2)                        // Z2Z2 = Z2^2 (mag: 1)
	u1.Set(x1).Mul(&z2z2).Normalize()         // U1 = X1*Z2Z2 (mag: 1)
	u2.Set(x2).Mul(&z1z1).Normalize()         // U2 = X2*Z1Z1 (mag: 1)
	s1.Set(y1).Mul(&z2z2).Mul(z2).Normalize() // S1 = Y1*Z2*Z2Z2 (mag: 1)
	s2.Set(y2).Mul(&z1z1).Mul(z1).Normalize() // S2 = Y2*Z1*Z1Z1 (mag: 1)
	if u1.Equals(&u2) {
		if s1.Equals(&s2) {
			// Since x1 == x2 and y1 == y2, point doubling must be
			// done, otherwise the addition would end up dividing
			// by zero.
			curve.doubleJacobian(x1, y1, z1, x3, y3, z3)
			return
		}

		// Since x1 == x2 and y1 == -y2, the sum is the point at
		// infinity per the group law.
		x3.SetInt(0)
		y3.SetInt(0)
		z3.SetInt(0)
		return
	}

	// Calculate X3, Y3, and Z3 according to the intermediate elements
	// breakdown above.
	var h, i, j, r, rr, v FieldVal
	var negU1, negS1, negX3 FieldVal
	negU1.Set(&u1).Negate(1)               // negU1 = -U1 (mag: 2)
	h.Add2(&u2, &negU1)                    // H = U2-U1 (mag: 3)
	i.Set(&h).MulInt(2).Square()           // I = (2*H)^2 (mag: 2)
	j.Mul2(&h, &i)                         // J = H*I (mag: 1)
	negS1.Set(&s1).Negate(1)               // negS1 = -S1 (mag: 2)
	r.Set(&s2).Add(&negS1).MulInt(2)       // r = 2*(S2-S1) (mag: 6)
	rr.SquareVal(&r)                       // rr = r^2 (mag: 1)
	v.Mul2(&u1, &i)                        // V = U1*I (mag: 1)
	x3.Set(&v).MulInt(2).Add(&j).Negate(3) // X3 = -(J+2*V) (mag: 4)
	x3.Add(&rr)                            // X3 = r^2+X3 (mag: 5)
	negX3.Set(x3).Negate(5)                // negX3 = -X3 (mag: 6)
	y3.Mul2(&s1, &j).MulInt(2).Negate(2)   // Y3 = -(2*S1*J) (mag: 3)
	y3.Add(v.Add(&negX3).Mul(&r))          // Y3 = r*(V-X3)+Y3 (mag: 4)
	z3.Add2(z1, z2).Square()               // Z3 = (Z1+Z2)^2 (mag: 1)
	z3.Add(z1z1.Add(&z2z2).Negate(2))      // Z3 = Z3-(Z1Z1+Z2Z2) (mag: 4)
	z3.Mul(&h)                             // Z3 = Z3*H (mag: 1)

	// Normalize the resulting field values to a magnitude of 1 as needed.
	x3.Normalize()
	y3.Normalize()
}

// 批量加法，传入参数为批量的两个需要加法的值，返回进行批量加法后的值 这里要进行并行加速
// func (curve *KoblitzCurve) addAffineBatch(x1, y1, x2, y2 []*big.Int) ([]*big.Int, []*big.Int) {
// 	//加法默认就是不同两个点的加法
// 	// A point at infinity is the identity according to the group law for
// 	// elliptic curve cryptography.  Thus, ∞ + P = P and P + ∞ = P.
// 	// if y1[0] == nil {
// 	// 	return x1, y1
// 	// }
// 	//定义结果变量
// 	resultX := make([]*big.Int, len(x1))
// 	resultY := make([]*big.Int, len(y1))
// 	k := make([]*big.Int, len(x1))
// 	//进行批量处理,对k x2 - x1的逆进行一个计算，计算完后才能计算后面的
// 	for i := 0; i < len(x1); i++ {
// 		// //首先对出现的单位元进行处理
// 		// if x1[i].Sign() == 0 && y1[i].Sign() == 0 {
// 		// 	resultX[i] = x2[i]
// 		// 	resultY[i] = y2[i]
// 		// 	k[i] = big.NewInt(1)
// 		// 	continue
// 		// }
// 		// if x2[i].Sign() == 0 && y2[i].Sign() == 0 {
// 		// 	resultX[i] = x1[i]
// 		// 	resultY[i] = y1[i]
// 		// 	k[i] = big.NewInt(1)
// 		// 	continue
// 		// }
// 		// A point at infinity is the identity according to the group law for
// 		// elliptic curve cryptography.  Thus, ∞ + P = P and P + ∞ = P.
// 		if (x1[i] == nil && y1[i] == nil) || (x1[i].Cmp(big.NewInt(0)) == 0 && y1[i].Cmp(big.NewInt(0)) == 0) {
// 			resultX[i] = x2[i]
// 			resultY[i] = y2[i]
// 			k[i] = big.NewInt(1)
// 			continue
// 		}
// 		if (x2[i] == nil && y2[i] == nil) || (x2[i].Cmp(big.NewInt(0)) == 0 && y2[i].Cmp(big.NewInt(0)) == 0) {
// 			resultX[i] = x1[i]
// 			resultY[i] = y1[i]
// 			k[i] = big.NewInt(1)
// 			continue
// 		}
// 		//计算x2-x1的值
// 		k[i] = new(big.Int).Sub(x2[i], x1[i])
// 	}

// 	//对k[i]中的元素进行Montgomery Trick的应用 进行求逆
// 	runtime.GOMAXPROCS(Threadnum)
// 	kTree := kBuildTree(k) //构建树 k数组就是x2-x1
// 	//treeroot_inv := new(big.Int).Set(kTree[Treelen-2]).Inverse() //获取树的根节点，然后计算其逆元
// 	treeroot_inv := new(big.Int).ModInverse(kTree[Treelen-2], c.P) //获取树的根节点，然后计算其逆元
// 	ZinvTree := kGetInvTree(treeroot_inv, kTree)                   //构建逆树

// 	//并行加速单个加法阶段
// 	overthreadnum := make(chan int, Threadnum) //创建和线程数相同的缓冲通道
// 	batch := Imax / (Threadnum)                //分组数

// 	for t := 0; t < Threadnum; t++ { //t循环迭代的编号， batch循环迭代的范围
// 		//利用逆树进行加法的函数
// 		go Get3(x1, y1, x2, y2, resultX, resultY, ZinvTree, t*batch, (t+1)*batch, overthreadnum) //利用Montgomery Trick 来批量加法计算
// 	}

// 	for i := 0; i < Threadnum; i++ {
// 		<-overthreadnum
// 	}
// 	WG.Done()

// 	return resultX, resultY

// }

// 不要并行
func (curve *KoblitzCurve) addAffineBatch(x1, y1, x2, y2 []*big.Int) ([]*big.Int, []*big.Int) {
	//加法默认就是不同两个点的加法
	// A point at infinity is the identity according to the group law for
	// elliptic curve cryptography.  Thus, ∞ + P = P and P + ∞ = P.
	// if y1[0] == nil {
	// 	return x1, y1
	// }
	//定义结果变量
	resultX := make([]*big.Int, len(x1))
	resultY := make([]*big.Int, len(y1))
	// 初始化 resultX 和 resultY 中的每个元素为一个新的 big.Int 变量
	for i := range resultX {
		resultX[i] = new(big.Int)
	}

	for i := range resultY {
		resultY[i] = new(big.Int)
	}
	k := make([]*big.Int, len(x1))
	//进行批量处理,对k x2 - x1的逆进行一个计算，计算完后才能计算后面的
	for i := 0; i < len(x1); i++ {
		// //首先对出现的单位元进行处理
		// if x1[i].Sign() == 0 && y1[i].Sign() == 0 {
		// 	resultX[i] = x2[i]
		// 	resultY[i] = y2[i]
		// 	k[i] = big.NewInt(1)
		// 	continue
		// }
		// if x2[i].Sign() == 0 && y2[i].Sign() == 0 {
		// 	resultX[i] = x1[i]
		// 	resultY[i] = y1[i]
		// 	k[i] = big.NewInt(1)
		// 	continue
		// }
		// A point at infinity is the identity according to the group law for
		// elliptic curve cryptography.  Thus, ∞ + P = P and P + ∞ = P.
		if (x1[i] == nil && y1[i] == nil) || (x1[i].Cmp(big.NewInt(0)) == 0 && y1[i].Cmp(big.NewInt(0)) == 0) {
			resultX[i] = x2[i]
			resultY[i] = y2[i]
			k[i] = big.NewInt(1)
			continue
		}
		if (x2[i] == nil && y2[i] == nil) || (x2[i].Cmp(big.NewInt(0)) == 0 && y2[i].Cmp(big.NewInt(0)) == 0) {
			resultX[i] = x1[i]
			resultY[i] = y1[i]
			k[i] = big.NewInt(1)
			continue
		}
		//计算x2-x1的值
		k[i] = new(big.Int).Sub(x2[i], x1[i])
	}

	//对k[i]中的元素进行Montgomery Trick的应用 进行求逆

	kTree := kBuildTree(k) //构建树 k数组就是x2-x1
	//treeroot_inv := new(big.Int).Set(kTree[Treelen-2]).Inverse() //获取树的根节点，然后计算其逆元
	treeroot_inv := new(big.Int).ModInverse(kTree[Treelen-2], c.P) //获取树的根节点，然后计算其逆元
	ZinvTree := kGetInvTree(treeroot_inv, kTree)                   //构建逆树

	//并行加速单个加法阶段
	overthreadnum := make(chan int, Threadnum) //创建和线程数相同的缓冲通道

	//利用逆树进行加法的函数
	Get3(x1, y1, x2, y2, resultX, resultY, ZinvTree, 0, Imax, overthreadnum) //利用Montgomery Trick 来批量加法计算

	return resultX, resultY

}

// 不要并行 FieldVal 类型上
func (curve *KoblitzCurve) faddAffineBatch(x1, y1, x2, y2, x3, y3 []*FieldVal) {
	//加法默认就是不同两个点的加法
	// A point at infinity is the identity according to the group law for
	// elliptic curve cryptography.  Thus, ∞ + P = P and P + ∞ = P.
	// if y1[0] == nil {
	// 	return x1, y1
	// }
	var k []*FieldVal
	//var s FieldVal
	//k := make([]*FieldVal, len(x1))
	// 初始化 x2, y2 中的每个元素为一个新的 FieldVal 变量
	// for i := range x1 {
	// 	k[i] = new(FieldVal)
	// }

	//进行批量处理,对k x2 - x1的逆进行一个计算，计算完后才能计算后面的
	for i := 0; i < len(x1); i++ {
		k = append(k, new(FieldVal))
		//k[i] = new(FieldVal)
		// //首先对出现的单位元进行处理
		// if x1[i].Sign() == 0 && y1[i].Sign() == 0 {
		// 	resultX[i] = x2[i]
		// 	resultY[i] = y2[i]
		// 	k[i] = big.NewInt(1)
		// 	continue
		// }
		// if x2[i].Sign() == 0 && y2[i].Sign() == 0 {
		// 	resultX[i] = x1[i]
		// 	resultY[i] = y1[i]
		// 	k[i] = big.NewInt(1)
		// 	continue
		// }
		// A point at infinity is the identity according to the group law for
		// elliptic curve cryptography.  Thus, ∞ + P = P and P + ∞ = P.
		if x1[i].IsZero() && y1[i].IsZero() {
			x3[i].Set(x2[i])
			y3[i].Set(y2[i])
			k[i].SetInt(1)
			continue
		} else if x2[i].IsZero() && y2[i].IsZero() {
			x3[i].Set(x1[i])
			y3[i].Set(y1[i])
			k[i].SetInt(1)
			continue
		}
		//fmt.Println("--------------------------------")

		//计算x2-x1的值
		var s FieldVal
		//fmt.Println(x1[i], x2[i])
		s.Set(x1[i]).Negate(1)
		//fmt.Println("s", s)
		//s.Normalize()
		//fmt.Println("s.Normalize()", s)
		k[i].Set(&s).Add(x2[i])
		//fmt.Println("k[i].Set(&s).Add(x2[i])", k[i])
		//k[i].Normalize()
		//fmt.Println("k[i].Normalize()", k[i])
		//k[i].Set((x1[i]).Negate(1)).Add(x2[i])

	}
	// fmt.Println("--------------------------------")
	// fmt.Println(x1[0], x2[0])
	// fmt.Println(k[0])
	// fmt.Println(k[0], "k[0]  x2-x1的值----1")

	//fmt.Println((k[0]).Inverse(), "k[0]  x2-x1逆的值----2") //这个语句会把k[0]元素的值取逆再赋值
	//对k[i]中的元素进行Montgomery Trick的应用 进行求逆

	//start3 := time.Now()
	kTree := BuildTree(k) //构建树 k数组就是x2-x1
	// fmt.Println(kTree[0], "kTree[0]  x2 -x1的值----2")
	// fmt.Println((kTree[0]).Inverse(), "kTree[0]  x2 -x1逆的值----3")
	//treeroot_inv := new(big.Int).Set(kTree[Treelen-2]).Inverse() //获取树的根节点，然后计算其逆元
	treeroot_inv := new(FieldVal).Set(kTree[Treelen-2].Inverse()) //获取树的根节点，然后计算其逆元
	ZinvTree := GetInvTree(treeroot_inv, kTree)                   //构建逆树
	//cost3 := time.Since(start3)
	// fmt.Println(ZinvTree[0], "ZinvTree[0]  x2 -x1的逆值----3")
	// fmt.Println("--------------------------------")
	//逆树没有问题test1 := k[0].Inverse()
	// test1 := kTree[0].Inverse()
	// test2 := kTree[1].Inverse()
	// test3 := kTree[2].Inverse()
	// test4 := kTree[3].Inverse()
	// test5 := kTree[15].Inverse()

	// fmt.Println("--------------------------------")
	// fmt.Println(test1)
	// fmt.Println(ZinvTree[0])
	// fmt.Println(test2)
	// fmt.Println(ZinvTree[1])
	// fmt.Println(test3)
	// fmt.Println(ZinvTree[2])
	// fmt.Println(test4)
	// fmt.Println(ZinvTree[3])
	// fmt.Println(test5)
	// fmt.Println(ZinvTree[15])
	// fmt.Println("--------------------------------")
	//并行加速单个加法阶段
	overthreadnum := make(chan int, Threadnum) //创建和线程数相同的缓冲通道

	//利用逆树进行加法的函数
	//start4 := time.Now()
	fGet3(x1, y1, x2, y2, x3, y3, ZinvTree, 0, Imax, overthreadnum) //利用Montgomery Trick 来批量加法计算
	//cost4 := time.Since(start4)

	// fmt.Printf("kTree  cost=[%s]\n", cost3)
	// fmt.Printf("fGet2  cost=[%s]\n", cost4)
	// fmt.Printf("faddAffineBatch kTree %f times  fGet2\n", float64(cost3)/(float64(cost4)+float64(cost3)))
	// fmt.Printf("------------------------------\n")

}

// 不要并行 FieldVal 类型上
func (curve *KoblitzCurve) FaddAffineBatch(x1, y1, x2, y2, x3, y3 []*FieldVal) {
	//加法默认就是不同两个点的加法
	// A point at infinity is the identity according to the group law for
	// elliptic curve cryptography.  Thus, ∞ + P = P and P + ∞ = P.
	// if y1[0] == nil {
	// 	return x1, y1
	// }

	k := make([]*FieldVal, len(x1))
	// 初始化 x2, y2 中的每个元素为一个新的 FieldVal 变量
	// for i := range x1 {
	// 	k[i] = new(FieldVal)
	// }

	//进行批量处理,对k x2 - x1的逆进行一个计算，计算完后才能计算后面的
	for i := 0; i < len(x1); i++ {
		k[i] = new(FieldVal)
		// //首先对出现的单位元进行处理
		// if x1[i].Sign() == 0 && y1[i].Sign() == 0 {
		// 	resultX[i] = x2[i]
		// 	resultY[i] = y2[i]
		// 	k[i] = big.NewInt(1)
		// 	continue
		// }
		// if x2[i].Sign() == 0 && y2[i].Sign() == 0 {
		// 	resultX[i] = x1[i]
		// 	resultY[i] = y1[i]
		// 	k[i] = big.NewInt(1)
		// 	continue
		// }
		// A point at infinity is the identity according to the group law for
		// elliptic curve cryptography.  Thus, ∞ + P = P and P + ∞ = P.
		if x1[i].IsZero() && y1[i].IsZero() {
			x3[i].Set(x2[i])
			y3[i].Set(y2[i])
			k[i].SetInt(1)
			continue
		} else if x2[i].IsZero() && y2[i].IsZero() {
			x3[i].Set(x1[i])
			y3[i].Set(y1[i])
			k[i].SetInt(1)
			continue
		}
		//计算x2-x1的值
		k[i].Set(x1[i].Negate(1)).Add(x2[i])

	}

	//对k[i]中的元素进行Montgomery Trick的应用 进行求逆

	kTree := BuildTree(k) //构建树 k数组就是x2-x1
	//treeroot_inv := new(big.Int).Set(kTree[Treelen-2]).Inverse() //获取树的根节点，然后计算其逆元
	treeroot_inv := new(FieldVal).Set(kTree[Treelen-2].Inverse()) //获取树的根节点，然后计算其逆元
	ZinvTree := GetInvTree(treeroot_inv, kTree)                   //构建逆树

	//并行加速单个加法阶段
	overthreadnum := make(chan int, Threadnum) //创建和线程数相同的缓冲通道

	//利用逆树进行加法的函数
	fGet3(x1, y1, x2, y2, x3, y3, ZinvTree, 0, Imax, overthreadnum) //利用Montgomery Trick 来批量加法计算

}

// // 批量乘，传入参数为批量的一个需要倍乘的点坐标数组，返回进行批量倍乘后的点坐标数组 这里要进行并行加速
// func (curve *KoblitzCurve) doubleAffineBatch(x1, y1 []*big.Int) ([]*big.Int, []*big.Int) {

// 	//定义结果变量
// 	resultX := make([]*big.Int, len(x1))
// 	resultY := make([]*big.Int, len(y1))
// 	k := make([]*big.Int, len(x1))

// 	//计算倍乘中需要求逆的2y1，并且对所有点出现单位元或零元素的特殊情况进行处理排除
// 	for i := 0; i < len(x1); i++ {

// 		//if (x1[i] == nil && y1[i] == nil) || (x1[i].Cmp(big.NewInt(0)) == 0 && y1[i].Cmp(big.NewInt(0)) == 0) {
// 		if (x1[i] == nil && y1[i] == nil) || (x1[i].Cmp(big.NewInt(0)) == 0 && y1[i].Cmp(big.NewInt(0)) == 0) {
// 			resultX[i] = x1[i]
// 			resultY[i] = y1[i]
// 			k[i] = big.NewInt(1)
// 			continue
// 		}

// 		//k[i].Add(y1[i], y1[i]) // k[i]的值为2y1
// 		k[i] = new(big.Int).Add(y1[i], y1[i]) // k[i]的值为2y1
// 	}
// 	fmt.Println(k)
// 	fmt.Println(len(x1))

// 	//对k[i]中的元素进行Montgomery Trick的应用 进行求逆
// 	runtime.GOMAXPROCS(Threadnum)
// 	kTree := kBuildTree(k) //构建树 k数组就是2y1
// 	//treeroot_inv := new(big.Int).Set(kTree[Treelen-2]).Inverse() //获取树的根节点，然后计算其逆元
// 	treeroot_inv := new(big.Int).ModInverse(kTree[Treelen-2], c.P) //获取树的根节点，然后计算其逆元
// 	ZinvTree := kGetInvTree(treeroot_inv, kTree)                   //构建逆树

// 	//并行加速单个加法阶段
// 	overthreadnum := make(chan int, Threadnum) //创建和线程数相同的缓冲通道
// 	batch := Imax / (Threadnum)                //分组数

// 	fmt.Println(Imax)
// 	fmt.Println(Threadnum)
// 	fmt.Println(batch)
// 	for t := 0; t < Threadnum; t++ { //t循环迭代的编号， batch循环迭代的范围
// 		//利用逆树进行加法的函数
// 		go Get2(x1, y1, resultX, resultY, ZinvTree, t*batch, (t+1)*batch, overthreadnum) //利用Montgomery Trick 来批量加法计算
// 	}

// 	for i := 0; i < Threadnum; i++ {
// 		<-overthreadnum
// 	}
// 	WG.Done()

// 	return resultX, resultY

// }

// 不要并行
// 批量乘，传入参数为批量的一个需要倍乘的点坐标数组，返回进行批量倍乘后的点坐标数组 这里要进行并行加速
func (curve *KoblitzCurve) doubleAffineBatch(x1, y1 []*big.Int) ([]*big.Int, []*big.Int) {

	//定义结果变量
	resultX := make([]*big.Int, len(x1))
	resultY := make([]*big.Int, len(y1))

	// 初始化 resultX 和 resultY 中的每个元素为一个新的 big.Int 变量
	for i := range resultX {
		resultX[i] = new(big.Int)
	}

	for i := range resultY {
		resultY[i] = new(big.Int)
	}
	k := make([]*big.Int, len(x1))

	//计算倍乘中需要求逆的2y1，并且对所有点出现单位元或零元素的特殊情况进行处理排除
	for i := 0; i < len(x1); i++ {

		//if (x1[i] == nil && y1[i] == nil) || (x1[i].Cmp(big.NewInt(0)) == 0 && y1[i].Cmp(big.NewInt(0)) == 0) {
		if (x1[i] == nil && y1[i] == nil) || (x1[i].Cmp(big.NewInt(0)) == 0 && y1[i].Cmp(big.NewInt(0)) == 0) {
			resultX[i] = x1[i]
			resultY[i] = y1[i]
			k[i] = big.NewInt(1)
			continue
		}

		//k[i].Add(y1[i], y1[i]) // k[i]的值为2y1
		k[i] = new(big.Int).Add(y1[i], y1[i]) // k[i]的值为2y1
	}

	//对k[i]中的元素进行Montgomery Trick的应用 进行求逆

	kTree := kBuildTree(k) //构建树 k数组就是2y1
	//treeroot_inv := new(big.Int).Set(kTree[Treelen-2]).Inverse() //获取树的根节点，然后计算其逆元
	treeroot_inv := new(big.Int).ModInverse(kTree[Treelen-2], c.P) //获取树的根节点，然后计算其逆元
	ZinvTree := kGetInvTree(treeroot_inv, kTree)                   //构建逆树

	//并行加速单个加法阶段
	overthreadnum := make(chan int, Threadnum) //创建和线程数相同的缓冲通道
	//分组数

	// fmt.Println(Imax)
	// fmt.Println(Threadnum)

	//利用逆树进行加法的函数
	Get2(x1, y1, resultX, resultY, ZinvTree, 0, Imax, overthreadnum) //利用Montgomery Trick 来批量加法计算

	return resultX, resultY

}

// 不要并行 有限域版本
// 批量乘，传入参数为批量的一个需要倍乘的点坐标数组，返回进行批量倍乘后的点坐标数组 这里要进行并行加速
func (curve *KoblitzCurve) fdoubleAffineBatch(x1, y1, x2, y2 []*FieldVal) {
	var k []*FieldVal
	//k := make([]*FieldVal, len(x1))
	// // 初始化 x2, y2 中的每个元素为一个新的 FieldVal 变量
	// for i := range x1 {
	// 	k[i] = new(FieldVal)
	// }

	//计算倍乘中需要求逆的2y1，并且对所有点出现单位元或零元素的特殊情况进行处理排除
	for i := 0; i < len(x1); i++ {
		k = append(k, new(FieldVal))
		//k[i] = new(FieldVal)
		//if (x1[i] == nil && y1[i] == nil) || (x1[i].Cmp(big.NewInt(0)) == 0 && y1[i].Cmp(big.NewInt(0)) == 0) {
		if x1[i].IsZero() && y1[i].IsZero() {
			x2[i].Set(x1[i])
			y2[i].Set(y1[i])

			k[i].SetInt(1)
			continue
		}

		//k[i].Add(y1[i], y1[i]) // k[i]的值为2y1
		k[i].Set(y1[i]).MulInt(2) // k[i]的值为2y1
	}

	//对k[i]中的元素进行Montgomery Trick的应用 进行求逆

	// fmt.Println(new(big.Int).SetBytes(k[0].Bytes()[:]), "k[0]  2y1的值----1")

	// fmt.Println(new(big.Int).SetBytes(k[0].Inverse().Bytes()[:]), "k[0]  2y1逆的值----2")
	//start3 := time.Now()
	kTree := BuildTree(k) //构建树 k数组就是2y1
	//cost5 := time.Since(start3)
	//treeroot_inv := new(big.Int).Set(kTree[Treelen-2]).Inverse() //获取树的根节点，然后计算其逆元
	//start5 := time.Now()
	treeroot_inv := new(FieldVal).Set(kTree[Treelen-2].Inverse()) //获取树的根节点，然后计算其逆元
	//cost6 := time.Since(start5)
	//start6 := time.Now()
	ZinvTree := GetInvTree(treeroot_inv, kTree) //构建逆树
	//cost7 := time.Since(start6)
	//cost3 := time.Since(start3)
	//逆树没有问题test1 := k[0].Inverse()
	// test1 := k[0].Inverse()
	// test2 := k[1].Inverse()
	// test3 := k[2].Inverse()
	// test4 := k[3].Inverse()
	// test5 := k[15].Inverse()

	// fmt.Println("--------------------------------")
	// fmt.Println(test1)
	// fmt.Println(ZinvTree[0])
	// fmt.Println(test2)
	// fmt.Println(ZinvTree[1])
	// fmt.Println(test3)
	// fmt.Println(ZinvTree[2])
	// fmt.Println(test4)
	// fmt.Println(ZinvTree[3])
	// fmt.Println(test5)
	// fmt.Println(ZinvTree[15])
	// fmt.Println("--------------------------------")

	//并行加速单个加法阶段
	overthreadnum := make(chan int, Threadnum) //创建和线程数相同的缓冲通道
	//分组数

	// fmt.Println(Imax)
	// fmt.Println(Threadnum)

	//利用逆树进行倍乘的函数
	//start4 := time.Now()
	fGet2(x1, y1, x2, y2, ZinvTree, 0, Imax, overthreadnum) //利用Montgomery Trick 来批量加法计算
	//cost4 := time.Since(start4)

	// fmt.Printf("kTree  cost=[%s]\n", cost3)
	// fmt.Printf("fGet2  cost=[%s]\n", cost4)

	// fmt.Printf("BuildTree  cost=[%s]\n", cost5)
	// fmt.Printf("treeroot_inv  cost=[%s]\n", cost6)
	// fmt.Printf("GetInvTree  cost=[%s]\n", cost7)

	// fmt.Printf("fdoubleAffineBatch kTree %f times  fGet2\n", float64(cost3)/(float64(cost4)+float64(cost3)))
	// fmt.Printf("------------------------------\n")
	// fmt.Printf("BuildTree treeroot_inv   GetInvTree %f  %f   %f times  kTree\n", float64(cost5)/float64(cost3), float64(cost6)/float64(cost3), float64(cost7)/float64(cost3))

}

// 不要并行 有限域版本
// 批量乘，传入参数为批量的一个需要倍乘的点坐标数组，返回进行批量倍乘后的点坐标数组 这里要进行并行加速
func (curve *KoblitzCurve) FdoubleAffineBatch(x1, y1, x2, y2 []*FieldVal) {

	k := make([]*FieldVal, len(x1))
	// // 初始化 x2, y2 中的每个元素为一个新的 FieldVal 变量
	// for i := range x1 {
	// 	k[i] = new(FieldVal)
	// }

	//计算倍乘中需要求逆的2y1，并且对所有点出现单位元或零元素的特殊情况进行处理排除
	for i := 0; i < len(x1); i++ {
		k[i] = new(FieldVal)
		//if (x1[i] == nil && y1[i] == nil) || (x1[i].Cmp(big.NewInt(0)) == 0 && y1[i].Cmp(big.NewInt(0)) == 0) {
		if x1[i].IsZero() && y1[i].IsZero() {
			x2[i].Set(x1[i])
			y2[i].Set(y1[i])

			k[i].SetInt(1)
			continue
		}

		//k[i].Add(y1[i], y1[i]) // k[i]的值为2y1
		k[i].Set(y1[i]).MulInt(2) // k[i]的值为2y1
	}

	//对k[i]中的元素进行Montgomery Trick的应用 进行求逆

	kTree := BuildTree(k) //构建树 k数组就是2y1
	//treeroot_inv := new(big.Int).Set(kTree[Treelen-2]).Inverse() //获取树的根节点，然后计算其逆元
	treeroot_inv := new(FieldVal).Set(kTree[Treelen-2].Inverse()) //获取树的根节点，然后计算其逆元
	ZinvTree := GetInvTree(treeroot_inv, kTree)                   //构建逆树

	//并行加速单个加法阶段
	overthreadnum := make(chan int, Threadnum) //创建和线程数相同的缓冲通道
	//分组数

	// fmt.Println(Imax)
	// fmt.Println(Threadnum)

	//利用逆树进行倍乘的函数
	fGet2(x1, y1, x2, y2, ZinvTree, 0, Imax, overthreadnum) //利用Montgomery Trick 来批量加法计算

}

// addJacobian adds the passed Jacobian points (x1, y1, z1) and (x2, y2, z2)
// together and stores the result in (x3, y3, z3).
// 雅可比坐标相加
func (curve *KoblitzCurve) addJacobian(x1, y1, z1, x2, y2, z2, x3, y3, z3 *FieldVal) {
	// A point at infinity is the identity according to the group law for
	// elliptic curve cryptography.  Thus, ∞ + P = P and P + ∞ = P.
	if (x1.IsZero() && y1.IsZero()) || z1.IsZero() {
		x3.Set(x2)
		y3.Set(y2)
		z3.Set(z2)
		return
	}
	if (x2.IsZero() && y2.IsZero()) || z2.IsZero() {
		x3.Set(x1)
		y3.Set(y1)
		z3.Set(z1)
		return
	}

	// Faster point addition can be achieved when certain assumptions are
	// met.  For example, when both points have the same z value, arithmetic
	// on the z values can be avoided.  This section thus checks for these
	// conditions and calls an appropriate add function which is accelerated
	// by using those assumptions.
	z1.Normalize()
	z2.Normalize()
	isZ1One := z1.Equals(fieldOne)
	isZ2One := z2.Equals(fieldOne)
	switch {
	case isZ1One && isZ2One:
		curve.addZ1AndZ2EqualsOne(x1, y1, z1, x2, y2, x3, y3, z3)
		return
	case z1.Equals(z2):
		curve.addZ1EqualsZ2(x1, y1, z1, x2, y2, x3, y3, z3)
		return
	case isZ2One:
		curve.addZ2EqualsOne(x1, y1, z1, x2, y2, x3, y3, z3)
		return
	}

	// None of the above assumptions are true, so fall back to generic
	// point addition.
	curve.addGeneric(x1, y1, z1, x2, y2, z2, x3, y3, z3)
}

// Add returns the sum of (x1,y1) and (x2,y2). Part of the elliptic.Curve
// interface.
// 两个访射坐标相加，本质是变成雅可比再变回来
func (curve *KoblitzCurve) Add(x1, y1, x2, y2 *big.Int) (*big.Int, *big.Int) {
	// A point at infinity is the identity according to the group law for
	// elliptic curve cryptography.  Thus, ∞ + P = P and P + ∞ = P.
	if x1.Sign() == 0 && y1.Sign() == 0 {
		return x2, y2
	}
	if x2.Sign() == 0 && y2.Sign() == 0 {
		return x1, y1
	}

	// Convert the affine coordinates from big integers to field values
	// and do the point addition in Jacobian projective space.
	fx1, fy1 := curve.bigAffineToField(x1, y1)
	fx2, fy2 := curve.bigAffineToField(x2, y2)
	fx3, fy3, fz3 := new(FieldVal), new(FieldVal), new(FieldVal)
	fOne := new(FieldVal).SetInt(1)
	curve.addJacobian(fx1, fy1, fOne, fx2, fy2, fOne, fx3, fy3, fz3)

	// Convert the Jacobian coordinate field values back to affine big
	// integers.
	return curve.fieldJacobianToBigAffine(fx3, fy3, fz3)
}

func (curve *KoblitzCurve) AddToF(x1, y1, x2, y2 *big.Int) (*FieldVal, *FieldVal) {

	// Convert the affine coordinates from big integers to field values
	// and do the point addition in Jacobian projective space.
	fx1, fy1 := curve.bigAffineToField(x1, y1)
	fx2, fy2 := curve.bigAffineToField(x2, y2)
	fx3, fy3, fz3 := new(FieldVal), new(FieldVal), new(FieldVal)
	fOne := new(FieldVal).SetInt(1)
	curve.addJacobian(fx1, fy1, fOne, fx2, fy2, fOne, fx3, fy3, fz3)
	// Convert the Jacobian coordinate field values back to affine big
	// integers.
	bx3, by3 := curve.fieldJacobianToBigAffine(fx3, fy3, fz3)
	fbx3, fby3 := curve.bigAffineToField(bx3, by3)
	return fbx3, fby3
}

func (curve *KoblitzCurve) Add1(x1, y1, x2, y2 *big.Int) (*FieldVal, *FieldVal, *FieldVal) {
	// A point at infinity is the identity according to the group law for
	// elliptic curve cryptography.  Thus, ∞ + P = P and P + ∞ = P.

	// Convert the affine coordinates from big integers to field values
	// and do the point addition in Jacobian projective space.
	fx1, fy1 := curve.bigAffineToField(x1, y1)
	fx2, fy2 := curve.bigAffineToField(x2, y2)
	fx3, fy3, fz3 := new(FieldVal), new(FieldVal), new(FieldVal)
	fOne := new(FieldVal).SetInt(1)
	curve.addJacobian(fx1, fy1, fOne, fx2, fy2, fOne, fx3, fy3, fz3)
	// Convert the Jacobian coordinate field values back to affine big
	// integers.
	return fx3, fy3, fz3
}

// doubleZ1EqualsOne performs point doubling on the passed Jacobian point
// when the point is already known to have a z value of 1 and stores
// the result in (x3, y3, z3).  That is to say (x3, y3, z3) = 2*(x1, y1, 1).  It
// performs faster point doubling than the generic routine since less arithmetic
// is needed due to the ability to avoid multiplication by the z value.
// 雅可比坐标翻倍，当z =1
func (curve *KoblitzCurve) doubleZ1EqualsOne(x1, y1, x3, y3, z3 *FieldVal) {
	// This function uses the assumptions that z1 is 1, thus the point
	// doubling formulas reduce to:
	//
	// X3 = (3*X1^2)^2 - 8*X1*Y1^2
	// Y3 = (3*X1^2)*(4*X1*Y1^2 - X3) - 8*Y1^4
	// Z3 = 2*Y1
	//
	// To compute the above efficiently, this implementation splits the
	// equation into intermediate elements which are used to minimize the
	// number of field multiplications in favor of field squarings which
	// are roughly 35% faster than field multiplications with the current
	// implementation at the time this was written.
	//
	// This uses a slightly modified version of the method shown at:
	// http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-mdbl-2007-bl
	//
	// In particular it performs the calculations using the following:
	// A = X1^2, B = Y1^2, C = B^2, D = 2*((X1+B)^2-A-C)
	// E = 3*A, F = E^2, X3 = F-2*D, Y3 = E*(D-X3)-8*C
	// Z3 = 2*Y1
	//
	// This results in a cost of 1 field multiplication, 5 field squarings,
	// 6 field additions, and 5 integer multiplications.
	var a, b, c, d, e, f FieldVal
	z3.Set(y1).MulInt(2)                     // Z3 = 2*Y1 (mag: 2)
	a.SquareVal(x1)                          // A = X1^2 (mag: 1)
	b.SquareVal(y1)                          // B = Y1^2 (mag: 1)
	c.SquareVal(&b)                          // C = B^2 (mag: 1)
	b.Add(x1).Square()                       // B = (X1+B)^2 (mag: 1)
	d.Set(&a).Add(&c).Negate(2)              // D = -(A+C) (mag: 3)
	d.Add(&b).MulInt(2)                      // D = 2*(B+D)(mag: 8)
	e.Set(&a).MulInt(3)                      // E = 3*A (mag: 3)
	f.SquareVal(&e)                          // F = E^2 (mag: 1)
	x3.Set(&d).MulInt(2).Negate(16)          // X3 = -(2*D) (mag: 17)
	x3.Add(&f)                               // X3 = F+X3 (mag: 18)
	f.Set(x3).Negate(18).Add(&d).Normalize() // F = D-X3 (mag: 1)
	y3.Set(&c).MulInt(8).Negate(8)           // Y3 = -(8*C) (mag: 9)
	y3.Add(f.Mul(&e))                        // Y3 = E*F+Y3 (mag: 10)

	// Normalize the field values back to a magnitude of 1.
	x3.Normalize()
	y3.Normalize()
	z3.Normalize()
}

// doubleGeneric performs point doubling on the passed Jacobian point without
// any assumptions about the z value and stores the result in (x3, y3, z3).
// That is to say (x3, y3, z3) = 2*(x1, y1, z1).  It is the slowest of the point
// doubling routines due to requiring the most arithmetic.
// 雅可比坐标翻倍
func (curve *KoblitzCurve) doubleGeneric(x1, y1, z1, x3, y3, z3 *FieldVal) {
	// Point doubling formula for Jacobian coordinates for the secp256k1
	// curve:
	// X3 = (3*X1^2)^2 - 8*X1*Y1^2
	// Y3 = (3*X1^2)*(4*X1*Y1^2 - X3) - 8*Y1^4
	// Z3 = 2*Y1*Z1
	//
	// To compute the above efficiently, this implementation splits the
	// equation into intermediate elements which are used to minimize the
	// number of field multiplications in favor of field squarings which
	// are roughly 35% faster than field multiplications with the current
	// implementation at the time this was written.
	//
	// This uses a slightly modified version of the method shown at:
	// http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l
	//
	// In particular it performs the calculations using the following:
	// A = X1^2, B = Y1^2, C = B^2, D = 2*((X1+B)^2-A-C)
	// E = 3*A, F = E^2, X3 = F-2*D, Y3 = E*(D-X3)-8*C
	// Z3 = 2*Y1*Z1
	//
	// This results in a cost of 1 field multiplication, 5 field squarings,
	// 6 field additions, and 5 integer multiplications.
	var a, b, c, d, e, f FieldVal
	z3.Mul2(y1, z1).MulInt(2)                // Z3 = 2*Y1*Z1 (mag: 2)
	a.SquareVal(x1)                          // A = X1^2 (mag: 1)
	b.SquareVal(y1)                          // B = Y1^2 (mag: 1)
	c.SquareVal(&b)                          // C = B^2 (mag: 1)
	b.Add(x1).Square()                       // B = (X1+B)^2 (mag: 1)
	d.Set(&a).Add(&c).Negate(2)              // D = -(A+C) (mag: 3)
	d.Add(&b).MulInt(2)                      // D = 2*(B+D)(mag: 8)
	e.Set(&a).MulInt(3)                      // E = 3*A (mag: 3)
	f.SquareVal(&e)                          // F = E^2 (mag: 1)
	x3.Set(&d).MulInt(2).Negate(16)          // X3 = -(2*D) (mag: 17)
	x3.Add(&f)                               // X3 = F+X3 (mag: 18)
	f.Set(x3).Negate(18).Add(&d).Normalize() // F = D-X3 (mag: 1)
	y3.Set(&c).MulInt(8).Negate(8)           // Y3 = -(8*C) (mag: 9)
	y3.Add(f.Mul(&e))                        // Y3 = E*F+Y3 (mag: 10)

	// Normalize the field values back to a magnitude of 1.
	x3.Normalize()
	y3.Normalize()
	z3.Normalize()
}

// doubleJacobian doubles the passed Jacobian point (x1, y1, z1) and stores the
// result in (x3, y3, z3).
func (curve *KoblitzCurve) doubleJacobian(x1, y1, z1, x3, y3, z3 *FieldVal) {
	// Doubling a point at infinity is still infinity.
	if y1.IsZero() || z1.IsZero() {
		x3.SetInt(0)
		y3.SetInt(0)
		z3.SetInt(0)
		return
	}

	// Slightly faster point doubling can be achieved when the z value is 1
	// by avoiding the multiplication on the z value.  This section calls
	// a point doubling function which is accelerated by using that
	// assumption when possible.
	if z1.Normalize().Equals(fieldOne) {
		curve.doubleZ1EqualsOne(x1, y1, x3, y3, z3)
		return
	}

	// Fall back to generic point doubling which works with arbitrary z
	// values.
	curve.doubleGeneric(x1, y1, z1, x3, y3, z3)
}

// Double returns 2*(x1,y1). Part of the elliptic.Curve interface.
func (curve *KoblitzCurve) Double(x1, y1 *big.Int) (*big.Int, *big.Int) {
	if y1.Sign() == 0 {
		return new(big.Int), new(big.Int)
	}

	// Convert the affine coordinates from big integers to field values
	// and do the point doubling in Jacobian projective space.
	fx1, fy1 := curve.bigAffineToField(x1, y1)
	fx3, fy3, fz3 := new(FieldVal), new(FieldVal), new(FieldVal)
	fOne := new(FieldVal).SetInt(1)
	curve.doubleJacobian(fx1, fy1, fOne, fx3, fy3, fz3)

	// Convert the Jacobian coordinate field values back to affine big
	// integers.
	return curve.fieldJacobianToBigAffine(fx3, fy3, fz3)
}

// splitK returns a balanced length-two representation of k and their signs.
// This is algorithm 3.74 from [GECC].
//
// One thing of note about this algorithm is that no matter what c1 and c2 are,
// the final equation of k = k1 + k2 * lambda (mod n) will hold.  This is
// provable mathematically due to how a1/b1/a2/b2 are computed.
//
// c1 and c2 are chosen to minimize the max(k1,k2).
func (curve *KoblitzCurve) splitK(k []byte) ([]byte, []byte, int, int) {
	// All math here is done with big.Int, which is slow.
	// At some point, it might be useful to write something similar to
	// FieldVal but for N instead of P as the prime field if this ends up
	// being a bottleneck.
	bigIntK := new(big.Int)
	c1, c2 := new(big.Int), new(big.Int)
	tmp1, tmp2 := new(big.Int), new(big.Int)
	k1, k2 := new(big.Int), new(big.Int)

	bigIntK.SetBytes(k)
	// c1 = round(b2 * k / n) from step 4.
	// Rounding isn't really necessary and costs too much, hence skipped
	c1.Mul(curve.b2, bigIntK)
	c1.Div(c1, curve.N)
	// c2 = round(b1 * k / n) from step 4 (sign reversed to optimize one step)
	// Rounding isn't really necessary and costs too much, hence skipped
	c2.Mul(curve.b1, bigIntK)
	c2.Div(c2, curve.N)
	// k1 = k - c1 * a1 - c2 * a2 from step 5 (note c2's sign is reversed)
	tmp1.Mul(c1, curve.a1)
	tmp2.Mul(c2, curve.a2)
	k1.Sub(bigIntK, tmp1)
	k1.Add(k1, tmp2)
	// k2 = - c1 * b1 - c2 * b2 from step 5 (note c2's sign is reversed)
	tmp1.Mul(c1, curve.b1)
	tmp2.Mul(c2, curve.b2)
	k2.Sub(tmp2, tmp1)

	// Note Bytes() throws out the sign of k1 and k2. This matters
	// since k1 and/or k2 can be negative. Hence, we pass that
	// back separately.
	return k1.Bytes(), k2.Bytes(), k1.Sign(), k2.Sign()
}

// moduloReduce reduces k from more than 32 bytes to 32 bytes and under.  This
// is done by doing a simple modulo curve.N.  We can do this since G^N = 1 and
// thus any other valid point on the elliptic curve has the same order.
func (curve *KoblitzCurve) moduloReduce(k []byte) []byte {
	// Since the order of G is curve.N, we can use a much smaller number
	// by doing modulo curve.N
	if len(k) > curve.byteSize {
		// Reduce k by performing modulo curve.N.
		tmpK := new(big.Int).SetBytes(k)
		tmpK.Mod(tmpK, curve.N)
		return tmpK.Bytes()
	}

	return k
}

// NAF takes a positive integer k and returns the Non-Adjacent Form (NAF) as two
// byte slices.  The first is where 1s will be.  The second is where -1s will
// be.  NAF is convenient in that on average, only 1/3rd of its values are
// non-zero.  This is algorithm 3.30 from [GECC].
//
// Essentially, this makes it possible to minimize the number of operations
// since the resulting ints returned will be at least 50% 0s.
func NAF(k []byte) ([]byte, []byte) {
	// The essence of this algorithm is that whenever we have consecutive 1s
	// in the binary, we want to put a -1 in the lowest bit and get a bunch
	// of 0s up to the highest bit of consecutive 1s.  This is due to this
	// identity:
	// 2^n + 2^(n-1) + 2^(n-2) + ... + 2^(n-k) = 2^(n+1) - 2^(n-k)
	//
	// The algorithm thus may need to go 1 more bit than the length of the
	// bits we actually have, hence bits being 1 bit longer than was
	// necessary.  Since we need to know whether adding will cause a carry,
	// we go from right-to-left in this addition.
	var carry, curIsOne, nextIsOne bool
	// these default to zero
	retPos := make([]byte, len(k)+1)
	retNeg := make([]byte, len(k)+1)
	for i := len(k) - 1; i >= 0; i-- {
		curByte := k[i]
		for j := uint(0); j < 8; j++ {
			curIsOne = curByte&1 == 1
			if j == 7 {
				if i == 0 {
					nextIsOne = false
				} else {
					nextIsOne = k[i-1]&1 == 1
				}
			} else {
				nextIsOne = curByte&2 == 2
			}
			if carry {
				if curIsOne {
					// This bit is 1, so continue to carry
					// and don't need to do anything.
				} else {
					// We've hit a 0 after some number of
					// 1s.
					if nextIsOne {
						// Start carrying again since
						// a new sequence of 1s is
						// starting.
						retNeg[i+1] += 1 << j
					} else {
						// Stop carrying since 1s have
						// stopped.
						carry = false
						retPos[i+1] += 1 << j
					}
				}
			} else if curIsOne {
				if nextIsOne {
					// If this is the start of at least 2
					// consecutive 1s, set the current one
					// to -1 and start carrying.
					retNeg[i+1] += 1 << j
					carry = true
				} else {
					// This is a singleton, not consecutive
					// 1s.
					retPos[i+1] += 1 << j
				}
			}
			curByte >>= 1
		}
	}
	if carry {
		retPos[0] = 1
		return retPos, retNeg
	}
	return retPos[1:], retNeg[1:]
}

// 实现将 FieldVal 类型转换为 *big.Int 类型的方法
func (f *FieldVal) ToBigInt() *big.Int {
	// 定义一个切片来存储 uint32 类型的字段值
	var nums []uint32
	// 将 uint32 类型的字段值复制到切片中
	for _, v := range f.n {
		nums = append(nums, v)
	}
	// 将切片转换为 *big.Int 类型的值
	betaBigInt := new(big.Int).SetUint64(uint64(nums[0]) | uint64(nums[1])<<26 | uint64(nums[2])<<52 | uint64(nums[3])<<78 | uint64(nums[4])<<104 | uint64(nums[5])<<130 | uint64(nums[6])<<156 | uint64(nums[7])<<182 | uint64(nums[8])<<208 | uint64(nums[9])<<234)

	return betaBigInt
}

// 批量乘算法，采用仿射坐标的运算
func (curve *KoblitzCurve) ScalarMultBatch(Bx, By []*big.Int, k []byte) ([]*big.Int, []*big.Int) {
	// 批量传入对应的点坐标x，y和标量切片数组。返回的值为x,y坐标数组

	// 初始化结果切片
	qx := make([]*big.Int, len(Bx))
	qy := make([]*big.Int, len(By))

	p1x := make([]*big.Int, len(Bx))
	p1y := make([]*big.Int, len(By))
	p1yNeg := make([]*big.Int, len(By))
	p2x := make([]*big.Int, len(Bx))
	p2y := make([]*big.Int, len(By))
	p2yNeg := make([]*big.Int, len(By))

	//标量k进行分解为k1 和 k2 并记录他们的符号 signK1和signK2
	k1, k2, signK1, signK2 := curve.splitK(curve.moduloReduce(k))
	// fmt.Println("-----------------------")
	// fmt.Println(k)
	// fmt.Println(k1, k2, signK1, signK2)
	//进行批量处理,对k x2 - x1的逆进行一个计算，计算完后才能计算后面的

	// 将 curve.beta 转换为 *big.Int 类型的值
	betaBigInt := curve.beta.ToBigInt()
	for i := 0; i < len(Bx); i++ {
		//将Bx，By的值赋给p1x， p2y，并计算p1y的相反数p1yNeg
		p1x[i] = Bx[i]
		p1y[i] = By[i]
		p1yNeg[i] = new(big.Int).Neg(new(big.Int).Set(By[i]))

		// 计算 p2x，即 p1x 乘以 curve.beta
		// p2xf, _ := curve.bigAffineToField(p1x[i], p1y[i])
		// p2x[i] = new(FieldVal).Mul2(p2xf, curve.beta).ToBigInt()

		p2x[i] = new(big.Int).Mul(new(big.Int).Set(p1x[i]), betaBigInt)
		p2x[i].Mod(p2x[i], c.P)
		p2y[i] = p1y[i]
		p2yNeg[i] = new(big.Int).Neg(new(big.Int).Set(p2y[i]))
	}

	fmt.Println("1-----------------------") //p2[x]算出来不对  对应乘betaBigInt结果是错的。
	// fmt.Println(betaBigInt)
	// fmt.Println(curve.beta)
	// fmt.Println(p1x[0], p1y[0], p1yNeg[0])
	fmt.Println(p2x[0], p2y[0], p2yNeg[0])

	// fmt.Println(k)
	// fmt.Println(Bx, By)
	// fmt.Println(p1x, p1y, p1yNeg)
	if signK1 == -1 {
		p1y, p1yNeg = p1yNeg, p1y
	}
	if signK2 == -1 {
		p2y, p2yNeg = p2yNeg, p2y
	}

	//以标量乘算法的顺序，并行的是每一步同时计算批量个

	//NAF 表示是一种用于优化标量乘法运算的表示方式，在某些情况下，它可以减少乘法运算的次数。因此，NAF 表示通常会包含很多零。
	// NAF versions of k1 and k2 should have a lot more zeros.
	//
	// The Pos version of the bytes contain the +1s and the Neg versions
	// contain the -1s.
	k1PosNAF, k1NegNAF := NAF(k1)
	k2PosNAF, k2NegNAF := NAF(k2)
	k1Len := len(k1PosNAF)
	k2Len := len(k2PosNAF)

	m := k1Len
	if m < k2Len {
		m = k2Len
	}

	//从左到右的标量乘法 他通过遍历m次循环来执行标量乘法的迭代过程，其中m是k1和k2中非零位的最大数量
	// Add left-to-right using the NAF optimization.  See algorithm 3.77
	// from [GECC].  This should be faster overall since there will be a lot
	// more instances of 0, hence reducing the number of Jacobian additions
	// at the cost of 1 possible extra doubling.
	//fmt.Println(m)
	var k1BytePos, k1ByteNeg, k2BytePos, k2ByteNeg byte
	for i := 0; i < m; i++ {
		// Since we're going left-to-right, pad the front with 0s.
		//判断对k1和k2进行零比特填充
		if i < m-k1Len {
			k1BytePos = 0
			k1ByteNeg = 0
		} else {
			k1BytePos = k1PosNAF[i-(m-k1Len)]
			k1ByteNeg = k1NegNAF[i-(m-k1Len)]
		}
		if i < m-k2Len {
			k2BytePos = 0
			k2ByteNeg = 0
		} else {
			k2BytePos = k2PosNAF[i-(m-k2Len)]
			k2ByteNeg = k2NegNAF[i-(m-k2Len)]
		}
		//fmt.Println(i)
		//雅可比执行double-and-add算法
		for j := 7; j >= 0; j-- {
			// Q = 2 * Q 每一步直接执行倍乘
			// fmt.Println(qx)
			// fmt.Println(qy)
			qx, qy = curve.doubleAffineBatch(qx, qy)

			//根据每一位是否为1进行add操作，并且还要判断y的正负
			if k1BytePos&0x80 == 0x80 {
				qx, qy = curve.addAffineBatch(qx, qy, p1x, p1y)
			} else if k1ByteNeg&0x80 == 0x80 {
				qx, qy = curve.addAffineBatch(qx, qy, p1x, p1yNeg)
			}

			if k2BytePos&0x80 == 0x80 {
				qx, qy = curve.addAffineBatch(qx, qy, p2x, p2y)
			} else if k2ByteNeg&0x80 == 0x80 {
				qx, qy = curve.addAffineBatch(qx, qy, p2x, p2yNeg)
			}
			k1BytePos <<= 1
			k1ByteNeg <<= 1
			k2BytePos <<= 1
			k2ByteNeg <<= 1
		}
	}

	return qx, qy

}

func (curve *KoblitzCurve) FScalarMultBatch(Bx, By []*big.Int, k []byte) ([]*FieldVal, []*FieldVal) {
	// 批量传入对应的点坐标x，y和标量切片数组。返回的值为x,y坐标数组

	// fmt.Printf("other  cost=[%s]\n", cost1)
	// fmt.Printf("fieldJacobianToBigAffine  cost=[%s]\n", cost2)
	// fmt.Printf("fieldJacobianToBigAffine %f times other\n", float64(cost2)/(float64(cost1)+float64(cost2)))
	// 初始Q点
	// 初始化 qx 和 qy 数组 以及计算需要用到的数据
	// var qx, qy, p1x, p1y, p1yNeg, p2x, p2y, p2yNeg []*FieldVal
	// length := len(Bx)
	// qx = make([]*FieldVal, length)
	// qy = make([]*FieldVal, length)
	// p1x = make([]*FieldVal, length)
	// p1y = make([]*FieldVal, length)
	// p1yNeg = make([]*FieldVal, length)
	// p2x = make([]*FieldVal, length)
	// p2y = make([]*FieldVal, length)
	// p2yNeg = make([]*FieldVal, length)

	// for i := range Bx {
	// 	qx[i] = new(FieldVal)
	// 	qy[i] = new(FieldVal)
	// 	p1x[i] = new(FieldVal)
	// 	p1y[i] = new(FieldVal)
	// 	p1yNeg[i] = new(FieldVal)
	// 	p2x[i] = new(FieldVal)
	// 	p2y[i] = new(FieldVal)
	// 	p2yNeg[i] = new(FieldVal)
	// }

	var qx, qy, p1x, p1y, p1yNeg, p2x, p2y, p2yNeg []*FieldVal

	// // 创建切片并初始化每个元素
	// for range Bx {
	// 	qx = append(qx, new(FieldVal))
	// 	qy = append(qy, new(FieldVal))
	// 	p1x = append(p1x, new(FieldVal))
	// 	p1y = append(p1y, new(FieldVal))
	// 	p1yNeg = append(p1yNeg, new(FieldVal))
	// 	p2x = append(p2x, new(FieldVal))
	// 	p2y = append(p2y, new(FieldVal))
	// 	p2yNeg = append(p2yNeg, new(FieldVal))
	// }
	//标量k进行分解为k1 和 k2 并记录他们的符号 signK1和signK2
	k1, k2, signK1, signK2 := curve.splitK(curve.moduloReduce(k))
	// fmt.Println("----------FScalarMultBatch-------------")
	// fmt.Println(Bx[0], By[0], k, k1, k2, signK1, signK2)
	// fmt.Println(k1, k2, signK1, signK2)
	//进行批量处理,对k x2 - x1的逆进行一个计算，计算完后才能计算后面的

	for i := 0; i < len(Bx); i++ {
		//创建切片并初始化每个元素
		qx = append(qx, new(FieldVal))
		qy = append(qy, new(FieldVal))
		p1x = append(p1x, new(FieldVal))
		p1y = append(p1y, new(FieldVal))
		p1yNeg = append(p1yNeg, new(FieldVal))
		p2x = append(p2x, new(FieldVal))
		p2y = append(p2y, new(FieldVal))
		p2yNeg = append(p2yNeg, new(FieldVal))

		//将Bx，By的值赋给p1x， p2y，并计算p1y的相反数p1yNeg
		p1x[i], p1y[i] = curve.bigAffineToField(Bx[i], By[i])
		p1yNeg[i].NegateVal(p1y[i], 1)

		p2x[i].Mul2(p1x[i], curve.beta)
		p2y[i].Set(p1y[i])
		p2yNeg[i].NegateVal(p2y[i], 1)

	}
	// fmt.Println("----------FScalarMultBatch-------------")
	// fmt.Println(p1x[0], p1y[0], p1yNeg[0], p2x[0], p2y[0], p2yNeg[0])

	// fmt.Println("1-----------------------") //p2[x]算出来不对  对应乘betaBigInt结果是错的。
	// // fmt.Println(betaBigInt)
	// // fmt.Println(curve.beta)
	// // fmt.Println(p1x[0], p1y[0], p1yNeg[0])
	// fmt.Println(p2x[0], p2y[0], p2yNeg[0])

	// fmt.Println(k)
	// fmt.Println(Bx, By)
	// fmt.Println(p1x, p1y, p1yNeg)
	if signK1 == -1 {
		p1y, p1yNeg = p1yNeg, p1y
	}
	if signK2 == -1 {
		p2y, p2yNeg = p2yNeg, p2y
	}

	//以标量乘算法的顺序，并行的是每一步同时计算批量个

	//NAF 表示是一种用于优化标量乘法运算的表示方式，在某些情况下，它可以减少乘法运算的次数。因此，NAF 表示通常会包含很多零。
	// NAF versions of k1 and k2 should have a lot more zeros.
	//
	// The Pos version of the bytes contain the +1s and the Neg versions
	// contain the -1s.
	k1PosNAF, k1NegNAF := NAF(k1)
	k2PosNAF, k2NegNAF := NAF(k2)
	k1Len := len(k1PosNAF)
	k2Len := len(k2PosNAF)

	m := k1Len
	if m < k2Len {
		m = k2Len
	}

	// fmt.Println(k1PosNAF, k1NegNAF, k2PosNAF, k2NegNAF, k1Len, k2Len, m)
	//从左到右的标量乘法 他通过遍历m次循环来执行标量乘法的迭代过程，其中m是k1和k2中非零位的最大数量
	// Add left-to-right using the NAF optimization.  See algorithm 3.77
	// from [GECC].  This should be faster overall since there will be a lot
	// more instances of 0, hence reducing the number of Jacobian additions
	// at the cost of 1 possible extra doubling.
	//fmt.Println(m)
	var k1BytePos, k1ByteNeg, k2BytePos, k2ByteNeg byte
	for i := 0; i < m; i++ {
		// Since we're going left-to-right, pad the front with 0s.
		//判断对k1和k2进行零比特填充
		if i < m-k1Len {
			k1BytePos = 0
			k1ByteNeg = 0
		} else {
			k1BytePos = k1PosNAF[i-(m-k1Len)]
			k1ByteNeg = k1NegNAF[i-(m-k1Len)]
		}
		if i < m-k2Len {
			k2BytePos = 0
			k2ByteNeg = 0
		} else {
			k2BytePos = k2PosNAF[i-(m-k2Len)]
			k2ByteNeg = k2NegNAF[i-(m-k2Len)]
		}
		//fmt.Println(i)
		//雅可比执行double-and-add算法

		// t2qx, t2qy := curve.FieldToBigAffine(qx[0], qy[0])
		// t2p1x, t2p1y := curve.FieldToBigAffine(p1x[0], p1y[0])
		// t2p2x, t2p2y := curve.FieldToBigAffine(p2x[0], p2y[0])
		// fmt.Println("t2qx------------t2qy")
		// fmt.Println(t2qx, t2qy)
		// fmt.Println("t2p1x------------t2p1y")
		// fmt.Println(t2p1x, t2p1y)
		// fmt.Println("t2p2x------------t2p2y")
		// fmt.Println(t2p2x, t2p2y)
		for j := 7; j >= 0; j-- {
			// Q = 2 * Q 每一步直接执行倍乘
			// fmt.Println(qx)
			// fmt.Println(qy)
			// fmt.Println("----------FScalarMultBatch-------------")
			// fmt.Println(qx[0], qy[0])
			// t2x, t2y := curve.FieldToBigAffine(qx[0], qy[0])
			// fmt.Println(t2x, t2y)
			curve.fdoubleAffineBatch(qx, qy, qx, qy)
			// fmt.Println(curve.FieldToBigAffine(qx[0], qy[0]))
			// t1x, t1y := curve.FieldToBigAffine(qx[0], qy[0])
			// fmt.Println(t1x, t1y)
			// fmt.Println(curve.FieldToBigAffine(p1x[0], p1y[0]))
			// fmt.Println(qx[0], qy[0])
			//根据每一位是否为1进行add操作，并且还要判断y的正负
			if k1BytePos&0x80 == 0x80 {
				curve.faddAffineBatch(qx, qy, p1x, p1y, qx, qy)
				// fmt.Println("---222-------FScalarMultBatch-------------")
				// fmt.Println(curve.FieldToBigAffine(qx[0], qy[0]))
			} else if k1ByteNeg&0x80 == 0x80 {
				curve.faddAffineBatch(qx, qy, p1x, p1yNeg, qx, qy)

			}

			if k2BytePos&0x80 == 0x80 {
				curve.faddAffineBatch(qx, qy, p2x, p2y, qx, qy)

			} else if k2ByteNeg&0x80 == 0x80 {
				curve.faddAffineBatch(qx, qy, p2x, p2yNeg, qx, qy)

			}
			k1BytePos <<= 1
			k1ByteNeg <<= 1
			k2BytePos <<= 1
			k2ByteNeg <<= 1
		}
	}
	return qx, qy

}

var cost1 time.Duration = 0
var cost2 time.Duration = 0

// ScalarMult returns k*(Bx, By) where k is a big endian integer.
// Part of the elliptic.Curve interface.
// 用的雅可比最后返回的是变回仿射坐标的值
func (curve *KoblitzCurve) ScalarMult(Bx, By *big.Int, k []byte) (*big.Int, *big.Int) {
	// Point Q = ∞ (point at infinity).
	//初始Q点
	//fmt.Println(k)
	//start3 := time.Now()
	qx, qy, qz := new(FieldVal), new(FieldVal), new(FieldVal)

	// Decompose K into k1 and k2 in order to halve the number of EC ops.
	// See Algorithm 3.74 in [GECC].

	//标量k进行分解为k1 和 k2 并记录他们的符号 signK1和signK2
	k1, k2, signK1, signK2 := curve.splitK(curve.moduloReduce(k))
	// fmt.Println("----------ScalarMult-------------")
	// fmt.Println(Bx, By, k, k1, k2, signK1, signK2)
	// fmt.Println(k)
	// fmt.Println(k1, k2, signK1, signK2)
	// The main equation here to remember is:
	//   k * P = k1 * P + k2 * ϕ(P)
	//
	// P1 below is P in the equation, P2 below is ϕ(P) in the equation
	p1x, p1y := curve.bigAffineToField(Bx, By) //将仿射坐标的x，y转换为有限域上的元素
	p1yNeg := new(FieldVal).NegateVal(p1y, 1)  //计算y的相反数
	p1z := new(FieldVal).SetInt(1)             //设置z点 初试1

	// NOTE: ϕ(x,y) = (βx,y).  The Jacobian z coordinate is 1, so this math
	// goes through.
	p2x := new(FieldVal).Mul2(p1x, curve.beta)
	p2y := new(FieldVal).Set(p1y)
	p2yNeg := new(FieldVal).NegateVal(p2y, 1)
	p2z := new(FieldVal).SetInt(1)

	// fmt.Println(curve.beta.ToBigInt())
	// fmt.Println("----------ScalarMult-------------")
	// fmt.Println(p1x, p1y, p1yNeg, p2x, p2y, p2yNeg)

	// Flip the positive and negative values of the points as needed
	// depending on the signs of k1 and k2.  As mentioned in the equation
	// above, each of k1 and k2 are multiplied by the respective point.
	// Since -k * P is the same thing as k * -P, and the group law for
	// elliptic curves states that P(x, y) = -P(x, -y), it's faster and
	// simplifies the code to just make the point negative.
	if signK1 == -1 {
		p1y, p1yNeg = p1yNeg, p1y
	}
	if signK2 == -1 {
		p2y, p2yNeg = p2yNeg, p2y
	}

	//NAF 表示是一种用于优化标量乘法运算的表示方式，在某些情况下，它可以减少乘法运算的次数。因此，NAF 表示通常会包含很多零。
	// NAF versions of k1 and k2 should have a lot more zeros.
	//
	// The Pos version of the bytes contain the +1s and the Neg versions
	// contain the -1s.
	k1PosNAF, k1NegNAF := NAF(k1)
	k2PosNAF, k2NegNAF := NAF(k2)
	k1Len := len(k1PosNAF)
	k2Len := len(k2PosNAF)

	m := k1Len
	if m < k2Len {
		m = k2Len
	}
	// fmt.Println(k1PosNAF, k1NegNAF, k2PosNAF, k2NegNAF, k1Len, k2Len, m)
	//从左到右的标量乘法 他通过遍历m次循环来执行标量乘法的迭代过程，其中m是k1和k2中非零位的最大数量
	// Add left-to-right using the NAF optimization.  See algorithm 3.77
	// from [GECC].  This should be faster overall since there will be a lot
	// more instances of 0, hence reducing the number of Jacobian additions
	// at the cost of 1 possible extra doubling.
	var k1BytePos, k1ByteNeg, k2BytePos, k2ByteNeg byte
	for i := 0; i < m; i++ {
		// Since we're going left-to-right, pad the front with 0s.
		//判断对k1和k2进行零比特填充
		if i < m-k1Len {
			k1BytePos = 0
			k1ByteNeg = 0
		} else {
			k1BytePos = k1PosNAF[i-(m-k1Len)]
			k1ByteNeg = k1NegNAF[i-(m-k1Len)]
		}
		if i < m-k2Len {
			k2BytePos = 0
			k2ByteNeg = 0
		} else {
			k2BytePos = k2PosNAF[i-(m-k2Len)]
			k2ByteNeg = k2NegNAF[i-(m-k2Len)]
		}

		// t1qx, t1qy := curve.fieldJacobianToBigAffine(qx, qy, qz)
		// t1p1x, t1p1y := curve.fieldJacobianToBigAffine(p1x, p1y, p1z)
		// t1p2x, t1p2y := curve.fieldJacobianToBigAffine(p2x, p2y, p2z)
		// fmt.Println("t1qx------------t1qy")
		// fmt.Println(t1qx, t1qy)
		// fmt.Println("t1p1x------------t1p1y")
		// fmt.Println(t1p1x, t1p1y)
		// fmt.Println("t1p2x------------t1p2y")
		// fmt.Println(t1p2x, t1p2y)
		//雅可比执行double-and-add算法
		for j := 7; j >= 0; j-- {
			// Q = 2 * Q
			// fmt.Println("----------ScalarMult-------------")
			// t1x, t1y := curve.fieldJacobianToBigAffine(qx, qy, qz)
			// fmt.Println(t1x, t1y)
			curve.doubleJacobian(qx, qy, qz, qx, qy, qz)
			// fmt.Println("----------ScalarMult-------------")
			// fmt.Println(curve.fieldJacobianToBigAffine(qx, qy, qz))
			// t2x, t2y := curve.fieldJacobianToBigAffine(qx, qy, qz)
			// fmt.Println(t2x, t2y)
			// fmt.Println(curve.fieldJacobianToBigAffine(p1x, p1y, p1z))
			if k1BytePos&0x80 == 0x80 {
				curve.addJacobian(qx, qy, qz, p1x, p1y, p1z,
					qx, qy, qz)
				// fmt.Println("----------ScalarMult-------------")
				// fmt.Println(curve.fieldJacobianToBigAffine(qx, qy, qz))
				// fmt.Println("----------ScalarMult-------------")
				// fmt.Println(qx, qy)
			} else if k1ByteNeg&0x80 == 0x80 {
				curve.addJacobian(qx, qy, qz, p1x, p1yNeg, p1z,
					qx, qy, qz)

			}

			if k2BytePos&0x80 == 0x80 {
				curve.addJacobian(qx, qy, qz, p2x, p2y, p2z,
					qx, qy, qz)

			} else if k2ByteNeg&0x80 == 0x80 {
				curve.addJacobian(qx, qy, qz, p2x, p2yNeg, p2z,
					qx, qy, qz)

			}
			k1BytePos <<= 1
			k1ByteNeg <<= 1
			k2BytePos <<= 1
			k2ByteNeg <<= 1
		}
	}
	//cost3 := time.Since(start3)
	// Convert the Jacobian coordinate field values back to affine big.Ints.
	//start4 := time.Now()
	//qx1, qy1 := curve.fieldJacobianToBigAffine(qx, qy, qz)
	//cost4 := time.Since(start4)

	//cost1 += cost3
	//cost2 += cost4
	// fmt.Printf("other  cost=[%s]\n", cost1)
	// fmt.Printf("fieldJacobianToBigAffine  cost=[%s]\n", cost2)
	// fmt.Printf("fieldJacobianToBigAffine %f times other\n", float64(cost2)/(float64(cost1)+float64(cost2)))
	//return qx1, qy1
	return curve.fieldJacobianToBigAffine(qx, qy, qz)

}

func (curve *KoblitzCurve) scalarMult(Bx, By *big.Int, k []byte) (*FieldVal, *FieldVal, *FieldVal) {
	//start1 := time.Now()
	// Point Q = ∞ (point at infinity).
	//初始Q点
	//start3 := time.Now()
	qx, qy, qz := new(FieldVal), new(FieldVal), new(FieldVal)

	// Decompose K into k1 and k2 in order to halve the number of EC ops.
	// See Algorithm 3.74 in [GECC].

	//标量k进行分解为k1 和 k2 并记录他们的符号 signK1和signK2
	k1, k2, signK1, signK2 := curve.splitK(curve.moduloReduce(k))
	// fmt.Println("----------ScalarMult-------------")
	// fmt.Println(Bx, By, k, k1, k2, signK1, signK2)
	// fmt.Println(k)
	// fmt.Println(k1, k2, signK1, signK2)
	// The main equation here to remember is:
	//   k * P = k1 * P + k2 * ϕ(P)
	//
	// P1 below is P in the equation, P2 below is ϕ(P) in the equation
	p1x, p1y := curve.bigAffineToField(Bx, By) //将仿射坐标的x，y转换为有限域上的元素
	p1yNeg := new(FieldVal).NegateVal(p1y, 1)  //计算y的相反数
	p1z := new(FieldVal).SetInt(1)             //设置z点 初试1

	// NOTE: ϕ(x,y) = (βx,y).  The Jacobian z coordinate is 1, so this math
	// goes through.
	p2x := new(FieldVal).Mul2(p1x, curve.beta)
	p2y := new(FieldVal).Set(p1y)
	p2yNeg := new(FieldVal).NegateVal(p2y, 1)
	p2z := new(FieldVal).SetInt(1)

	// fmt.Println(curve.beta.ToBigInt())
	// fmt.Println("----------ScalarMult-------------")
	// fmt.Println(p1x, p1y, p1yNeg, p2x, p2y, p2yNeg)

	// Flip the positive and negative values of the points as needed
	// depending on the signs of k1 and k2.  As mentioned in the equation
	// above, each of k1 and k2 are multiplied by the respective point.
	// Since -k * P is the same thing as k * -P, and the group law for
	// elliptic curves states that P(x, y) = -P(x, -y), it's faster and
	// simplifies the code to just make the point negative.
	if signK1 == -1 {
		p1y, p1yNeg = p1yNeg, p1y
	}
	if signK2 == -1 {
		p2y, p2yNeg = p2yNeg, p2y
	}

	//NAF 表示是一种用于优化标量乘法运算的表示方式，在某些情况下，它可以减少乘法运算的次数。因此，NAF 表示通常会包含很多零。
	// NAF versions of k1 and k2 should have a lot more zeros.
	//
	// The Pos version of the bytes contain the +1s and the Neg versions
	// contain the -1s.
	k1PosNAF, k1NegNAF := NAF(k1)
	k2PosNAF, k2NegNAF := NAF(k2)
	k1Len := len(k1PosNAF)
	k2Len := len(k2PosNAF)

	m := k1Len
	if m < k2Len {
		m = k2Len
	}
	// fmt.Println(k1PosNAF, k1NegNAF, k2PosNAF, k2NegNAF, k1Len, k2Len, m)
	//从左到右的标量乘法 他通过遍历m次循环来执行标量乘法的迭代过程，其中m是k1和k2中非零位的最大数量
	// Add left-to-right using the NAF optimization.  See algorithm 3.77
	// from [GECC].  This should be faster overall since there will be a lot
	// more instances of 0, hence reducing the number of Jacobian additions
	// at the cost of 1 possible extra doubling.
	var k1BytePos, k1ByteNeg, k2BytePos, k2ByteNeg byte
	for i := 0; i < m; i++ {
		// Since we're going left-to-right, pad the front with 0s.
		//判断对k1和k2进行零比特填充
		if i < m-k1Len {
			k1BytePos = 0
			k1ByteNeg = 0
		} else {
			k1BytePos = k1PosNAF[i-(m-k1Len)]
			k1ByteNeg = k1NegNAF[i-(m-k1Len)]
		}
		if i < m-k2Len {
			k2BytePos = 0
			k2ByteNeg = 0
		} else {
			k2BytePos = k2PosNAF[i-(m-k2Len)]
			k2ByteNeg = k2NegNAF[i-(m-k2Len)]
		}

		// t1qx, t1qy := curve.fieldJacobianToBigAffine(qx, qy, qz)
		// t1p1x, t1p1y := curve.fieldJacobianToBigAffine(p1x, p1y, p1z)
		// t1p2x, t1p2y := curve.fieldJacobianToBigAffine(p2x, p2y, p2z)
		// fmt.Println("t1qx------------t1qy")
		// fmt.Println(t1qx, t1qy)
		// fmt.Println("t1p1x------------t1p1y")
		// fmt.Println(t1p1x, t1p1y)
		// fmt.Println("t1p2x------------t1p2y")
		// fmt.Println(t1p2x, t1p2y)
		//雅可比执行double-and-add算法
		for j := 7; j >= 0; j-- {
			// Q = 2 * Q
			// fmt.Println("----------ScalarMult-------------")
			// t1x, t1y := curve.fieldJacobianToBigAffine(qx, qy, qz)
			// fmt.Println(t1x, t1y)
			curve.doubleJacobian(qx, qy, qz, qx, qy, qz)
			// fmt.Println("----------ScalarMult-------------")
			// fmt.Println(curve.fieldJacobianToBigAffine(qx, qy, qz))
			// t2x, t2y := curve.fieldJacobianToBigAffine(qx, qy, qz)
			// fmt.Println(t2x, t2y)
			// fmt.Println(curve.fieldJacobianToBigAffine(p1x, p1y, p1z))
			if k1BytePos&0x80 == 0x80 {
				curve.addJacobian(qx, qy, qz, p1x, p1y, p1z,
					qx, qy, qz)
				// fmt.Println("----------ScalarMult-------------")
				// fmt.Println(curve.fieldJacobianToBigAffine(qx, qy, qz))
				// fmt.Println("----------ScalarMult-------------")
				// fmt.Println(qx, qy)
			} else if k1ByteNeg&0x80 == 0x80 {
				curve.addJacobian(qx, qy, qz, p1x, p1yNeg, p1z,
					qx, qy, qz)

			}

			if k2BytePos&0x80 == 0x80 {
				curve.addJacobian(qx, qy, qz, p2x, p2y, p2z,
					qx, qy, qz)

			} else if k2ByteNeg&0x80 == 0x80 {
				curve.addJacobian(qx, qy, qz, p2x, p2yNeg, p2z,
					qx, qy, qz)

			}
			k1BytePos <<= 1
			k1ByteNeg <<= 1
			k2BytePos <<= 1
			k2ByteNeg <<= 1
		}
	}
	///cost3 := time.Since(start3)
	// Convert the Jacobian coordinate field values back to affine big.Ints.
	//start4 := time.Now()
	//qx1, qy1 := curve.fieldJacobianToBigAffine(qx, qy, qz)
	//cost4 := time.Since(start4)
	// fmt.Printf("other  cost=[%s]\n", cost1)
	// fmt.Printf("fieldJacobianToBigAffine  cost=[%s]\n", cost2)
	// fmt.Printf("fieldJacobianToBigAffine %f times other\n", float64(cost2)/float64(cost1))
	//cost1 += cost3
	//cost2 += cost4
	//return qx1, qy1
	return qx, qy, qz

}

// 批量获得两不同点相加
func (curve *KoblitzCurve) PPGet3(x1, y1, x2, y2, P, Zinv *big.Int) (x *big.Int, y *big.Int) {
	//time.Sleep(time.Duration(4) * time.Second)

	var h, i, r, v, d, f, g, x3, y3 big.Int

	//fmt.Println("--------PPGet3------------")
	r.Sub(y2, y1) //r = y2 - y1
	//fmt.Println("r.Sub(y2, y1)", &r)
	//r.Mod(&r, P)
	//fmt.Println("r.Mod(&r, P)", &r)
	h.Mul(&r, Zinv) //h = (y2 - y1)/(x2 - x1)
	//fmt.Println("h.Mul(&r, Zinv)", &h)
	//h.Mod(&h, P) //h = (y2 - y1)/(x2 - x1) mod p
	//fmt.Println("h.Mod(&h, P)", &h)
	i.Exp(&h, big.NewInt(2), nil)
	//fmt.Println("i.Sqrt(&h)", &i)
	//i.Mod(&i, P)
	//fmt.Println("i.Mod(&i, P)", &i)
	v.Add(x1, x2).Neg(&v) //v = -x1 - x2
	//fmt.Println("v.Add(x1, x2).Neg(&v)", &v)
	f.Set(&i).Add(&i, &v) // f = h^2 - x1 - x2
	//fmt.Println("f.Set(&i).Add(&i, &v)", &f)
	x3.Mod(&f, P) // x3 = h^2 - x1 - x2 mod p
	//fmt.Println("x3.Mod(&f, P)", &x3)
	d.Sub(x1, &x3) //d = x1 -x3
	//fmt.Println("d.Sub(x1, &x3))", &d)
	g.Mul(&d, &h).Sub(&g, y1) // g = h*(x1 -x3) -y1
	//fmt.Println("g.Mul(&d, &h).Sub(&g, y1)", &g)
	y3.Mod(&g, P) //y3 = h*(x1 -x3) -y1 mod p
	//fmt.Println("y3.Mod(&g, P)", &y3)
	//fmt.Println("--------PPGet3------------")

	return &x3, &y3

}

func (curve *KoblitzCurve) ScalarMultBatchFieldVal(Bx, By []*big.Int, k []byte) ([]*big.Int, []*big.Int) {

	var resx, resy []*big.Int
	var Fx, Fy, Fz []*FieldVal
	//start1 := time.Now()
	for i := 0; i < len(Bx); i++ {
		//调用标量乘函数,返回还没进行最后一步转化
		Fx = append(Fx, new(FieldVal))
		Fy = append(Fy, new(FieldVal))
		Fz = append(Fz, new(FieldVal))
		Fx[i], Fy[i], Fz[i] = curve.scalarMult(Bx[i], By[i], k)
	}
	//costBatch1 := time.Since(start1)

	//使用trick
	//start2 := time.Now()
	kTree := BuildTree(Fz) //构建树 k数组就是
	//costBatchkTree := time.Since(start2)
	//start3 := time.Now()
	treeroot_inv := new(FieldVal).Set(kTree[Treelen-2].Inverse()) //获取树的根节点，然后计算其逆元
	//costBatchtreeroot_inv := time.Since(start3)
	//start4 := time.Now()
	ZinvTree := GetInvTree(treeroot_inv, kTree) //构建逆树
	//costBatchGetInvTree := time.Since(start4)

	//用逆树进行雅可比坐标的还原
	resx, resy = curve.fieldJacobianToBigAffineBatch(Fx, Fy, ZinvTree)

	// costBatch2 := time.Since(start2)
	// fmt.Printf("*******************************************\n")
	// fmt.Printf("ScalarMultBatchFieldVal的除开转化的用时 cost=[%s]\n", costBatch1)
	// fmt.Printf("ScalarMultBatchFieldVal的转化的用时 cost=[%s]\n", costBatch2)
	// fmt.Printf("ScalarMultBatchFieldVal的转化的建树用时 cost=[%s]\n", costBatchkTree)
	// fmt.Printf("ScalarMultBatchFieldVal的转化的求逆用时 cost=[%s]\n", costBatchtreeroot_inv)
	// fmt.Printf("ScalarMultBatchFieldVal的转化的求逆树用时 cost=[%s]\n", costBatchGetInvTree)
	// fmt.Printf("ScalarMultBatchFieldVal的转化的建树用时占转化部分 %f\n", float64(costBatchkTree)/float64(costBatch2))
	// fmt.Printf("ScalarMultBatchFieldVal的转化的求逆用时占转化部分 %f\n", float64(costBatchtreeroot_inv)/float64(costBatch2))
	// fmt.Printf("ScalarMultBatchFieldVal的转化的求逆树用时占转化部分 %f\n", float64(costBatchGetInvTree)/float64(costBatch2))
	// fmt.Printf("ScalarMultBatchFieldVal的转化的用时占比 %f\n", float64(costBatch2)/float64(costBatch1+costBatch2))

	// fmt.Printf("ScalarMult的除开转化的用时 cost=[%s]\n", cost1)
	// fmt.Printf("ScalarMult的转化的用时 cost=[%s]\n", cost2)
	// fmt.Printf("ScalarMult的转化的用时占比 %f\n", float64(cost2)/float64(cost1+cost2))
	// fmt.Printf("\n")
	// fmt.Printf("转化的用时优化率 %f\n", float64(cost2-costBatch2)/float64(cost2))
	// fmt.Printf("总用时优化率 %f\n", float64(cost1+cost2-costBatch1-costBatch2)/float64(cost1+cost2))
	// fmt.Printf("*******************************************\n")
	return resx, resy

}

func (curve *KoblitzCurve) ScalarMultBatchFieldValSingleQ(Bx, By *big.Int, ks [][]byte) ([]*big.Int, []*big.Int) {

	var resx, resy []*big.Int
	var Fx, Fy, Fz []*FieldVal

	for i := 0; i < Imax; i++ {
		//调用标量乘函数,返回还没进行最后一步转化
		Fx = append(Fx, new(FieldVal))
		Fy = append(Fy, new(FieldVal))
		Fz = append(Fz, new(FieldVal))
		Fx[i], Fy[i], Fz[i] = curve.scalarMult(Bx, By, ks[i])
	}

	//使用trick
	kTree := BuildTree(Fz) //构建树 k数组就是

	treeroot_inv := new(FieldVal).Set(kTree[Treelen-2].Inverse()) //获取树的根节点，然后计算其逆元

	ZinvTree := GetInvTree(treeroot_inv, kTree) //构建逆树

	//用逆树进行雅可比坐标的还原
	resx, resy = curve.fieldJacobianToBigAffineBatch(Fx, Fy, ZinvTree)

	return resx, resy

}

func (curve *KoblitzCurve) ScalarBaseMultBatchFieldVal(k [][]byte) ([]*big.Int, []*big.Int) {

	var resx, resy []*big.Int
	var Fx, Fy, Fz []*FieldVal

	for i := 0; i < Imax; i++ {
		//调用标量乘函数,返回还没进行最后一步转化
		Fx = append(Fx, new(FieldVal))
		Fy = append(Fy, new(FieldVal))
		Fz = append(Fz, new(FieldVal))
		Fx[i], Fy[i], Fz[i] = curve.ScalarBaseMultField(k[i])
	}

	//使用trick
	kTree := BuildTree(Fz) //构建树 k数组就是

	treeroot_inv := new(FieldVal).Set(kTree[Treelen-2].Inverse()) //获取树的根节点，然后计算其逆元

	ZinvTree := GetInvTree(treeroot_inv, kTree) //构建逆树

	//用逆树进行雅可比坐标的还原
	resx, resy = curve.fieldJacobianToBigAffineBatch(Fx, Fy, ZinvTree)

	return resx, resy

}

// ScalarMult returns k*(Bx, By) where k is a big endian integer.
// Part of the elliptic.Curve interface.
func (curve *KoblitzCurve) ScalarMultField(Bx, By *big.Int, k []byte) (*FieldVal, *FieldVal, *FieldVal) {
	//start1 := time.Now()
	// Point Q = ∞ (point at infinity).
	qx, qy, qz := new(FieldVal), new(FieldVal), new(FieldVal)

	// Decompose K into k1 and k2 in order to halve the number of EC ops.
	// See Algorithm 3.74 in [GECC].
	k1, k2, signK1, signK2 := curve.splitK(curve.moduloReduce(k))

	// The main equation here to remember is:
	//   k * P = k1 * P + k2 * ϕ(P)
	//
	// P1 below is P in the equation, P2 below is ϕ(P) in the equation
	p1x, p1y := curve.bigAffineToField(Bx, By)
	p1yNeg := new(FieldVal).NegateVal(p1y, 1)
	p1z := new(FieldVal).SetInt(1)

	// NOTE: ϕ(x,y) = (βx,y).  The Jacobian z coordinate is 1, so this math
	// goes through.
	p2x := new(FieldVal).Mul2(p1x, curve.beta)
	p2y := new(FieldVal).Set(p1y)
	p2yNeg := new(FieldVal).NegateVal(p2y, 1)
	p2z := new(FieldVal).SetInt(1)

	// Flip the positive and negative values of the points as needed
	// depending on the signs of k1 and k2.  As mentioned in the equation
	// above, each of k1 and k2 are multiplied by the respective point.
	// Since -k * P is the same thing as k * -P, and the group law for
	// elliptic curves states that P(x, y) = -P(x, -y), it's faster and
	// simplifies the code to just make the point negative.
	if signK1 == -1 {
		p1y, p1yNeg = p1yNeg, p1y
	}
	if signK2 == -1 {
		p2y, p2yNeg = p2yNeg, p2y
	}

	// NAF versions of k1 and k2 should have a lot more zeros.
	//
	// The Pos version of the bytes contain the +1s and the Neg versions
	// contain the -1s.
	k1PosNAF, k1NegNAF := NAF(k1)
	k2PosNAF, k2NegNAF := NAF(k2)
	k1Len := len(k1PosNAF)
	k2Len := len(k2PosNAF)

	m := k1Len
	if m < k2Len {
		m = k2Len
	}

	// Add left-to-right using the NAF optimization.  See algorithm 3.77
	// from [GECC].  This should be faster overall since there will be a lot
	// more instances of 0, hence reducing the number of Jacobian additions
	// at the cost of 1 possible extra doubling.
	var k1BytePos, k1ByteNeg, k2BytePos, k2ByteNeg byte
	for i := 0; i < m; i++ {
		// Since we're going left-to-right, pad the front with 0s.
		if i < m-k1Len {
			k1BytePos = 0
			k1ByteNeg = 0
		} else {
			k1BytePos = k1PosNAF[i-(m-k1Len)]
			k1ByteNeg = k1NegNAF[i-(m-k1Len)]
		}
		if i < m-k2Len {
			k2BytePos = 0
			k2ByteNeg = 0
		} else {
			k2BytePos = k2PosNAF[i-(m-k2Len)]
			k2ByteNeg = k2NegNAF[i-(m-k2Len)]
		}

		for j := 7; j >= 0; j-- {
			// Q = 2 * Q
			curve.doubleJacobian(qx, qy, qz, qx, qy, qz)

			if k1BytePos&0x80 == 0x80 {
				curve.addJacobian(qx, qy, qz, p1x, p1y, p1z,
					qx, qy, qz)
			} else if k1ByteNeg&0x80 == 0x80 {
				curve.addJacobian(qx, qy, qz, p1x, p1yNeg, p1z,
					qx, qy, qz)
			}

			if k2BytePos&0x80 == 0x80 {
				curve.addJacobian(qx, qy, qz, p2x, p2y, p2z,
					qx, qy, qz)
			} else if k2ByteNeg&0x80 == 0x80 {
				curve.addJacobian(qx, qy, qz, p2x, p2yNeg, p2z,
					qx, qy, qz)
			}
			k1BytePos <<= 1
			k1ByteNeg <<= 1
			k2BytePos <<= 1
			k2ByteNeg <<= 1
		}
	}

	// Convert the Jacobian coordinate field values back to affine big.Ints.
	return qx, qy, qz
}

func (curve *KoblitzCurve) ScalarMultPoint(P *Point, k []byte) *Point {
	//start1 := time.Now()
	// Point Q = ∞ (point at infinity).
	qx, qy, qz := new(FieldVal), new(FieldVal), new(FieldVal)

	// Decompose K into k1 and k2 in order to halve the number of EC ops.
	// See Algorithm 3.74 in [GECC].
	k1, k2, signK1, signK2 := curve.splitK(curve.moduloReduce(k))

	// The main equation here to remember is:
	//   k * P = k1 * P + k2 * ϕ(P)
	//
	// P1 below is P in the equation, P2 below is ϕ(P) in the equation
	p1x, p1y := curve.bigAffineToField(P.X, P.Y)
	p1yNeg := new(FieldVal).NegateVal(p1y, 1)
	p1z := new(FieldVal).SetInt(1)

	// NOTE: ϕ(x,y) = (βx,y).  The Jacobian z coordinate is 1, so this math
	// goes through.
	p2x := new(FieldVal).Mul2(p1x, curve.beta)
	p2y := new(FieldVal).Set(p1y)
	p2yNeg := new(FieldVal).NegateVal(p2y, 1)
	p2z := new(FieldVal).SetInt(1)

	// Flip the positive and negative values of the points as needed
	// depending on the signs of k1 and k2.  As mentioned in the equation
	// above, each of k1 and k2 are multiplied by the respective point.
	// Since -k * P is the same thing as k * -P, and the group law for
	// elliptic curves states that P(x, y) = -P(x, -y), it's faster and
	// simplifies the code to just make the point negative.
	if signK1 == -1 {
		p1y, p1yNeg = p1yNeg, p1y
	}
	if signK2 == -1 {
		p2y, p2yNeg = p2yNeg, p2y
	}

	// NAF versions of k1 and k2 should have a lot more zeros.
	//
	// The Pos version of the bytes contain the +1s and the Neg versions
	// contain the -1s.
	k1PosNAF, k1NegNAF := NAF(k1)
	k2PosNAF, k2NegNAF := NAF(k2)
	k1Len := len(k1PosNAF)
	k2Len := len(k2PosNAF)

	m := k1Len
	if m < k2Len {
		m = k2Len
	}

	// Add left-to-right using the NAF optimization.  See algorithm 3.77
	// from [GECC].  This should be faster overall since there will be a lot
	// more instances of 0, hence reducing the number of Jacobian additions
	// at the cost of 1 possible extra doubling.
	var k1BytePos, k1ByteNeg, k2BytePos, k2ByteNeg byte
	for i := 0; i < m; i++ {
		// Since we're going left-to-right, pad the front with 0s.
		if i < m-k1Len {
			k1BytePos = 0
			k1ByteNeg = 0
		} else {
			k1BytePos = k1PosNAF[i-(m-k1Len)]
			k1ByteNeg = k1NegNAF[i-(m-k1Len)]
		}
		if i < m-k2Len {
			k2BytePos = 0
			k2ByteNeg = 0
		} else {
			k2BytePos = k2PosNAF[i-(m-k2Len)]
			k2ByteNeg = k2NegNAF[i-(m-k2Len)]
		}

		for j := 7; j >= 0; j-- {
			// Q = 2 * Q
			curve.doubleJacobian(qx, qy, qz, qx, qy, qz)

			if k1BytePos&0x80 == 0x80 {
				curve.addJacobian(qx, qy, qz, p1x, p1y, p1z,
					qx, qy, qz)
			} else if k1ByteNeg&0x80 == 0x80 {
				curve.addJacobian(qx, qy, qz, p1x, p1yNeg, p1z,
					qx, qy, qz)
			}

			if k2BytePos&0x80 == 0x80 {
				curve.addJacobian(qx, qy, qz, p2x, p2y, p2z,
					qx, qy, qz)
			} else if k2ByteNeg&0x80 == 0x80 {
				curve.addJacobian(qx, qy, qz, p2x, p2yNeg, p2z,
					qx, qy, qz)
			}
			k1BytePos <<= 1
			k1ByteNeg <<= 1
			k2BytePos <<= 1
			k2ByteNeg <<= 1
		}
	}

	// Convert the Jacobian coordinate field values back to affine big.Ints.
	px, py := curve.fieldJacobianToBigAffine(qx, qy, qz)
	return &Point{px, py}
}

// ScalarBaseMult returns k*G where G is the base point of the group and k is a
// big endian integer.
// Part of the elliptic.Curve interface.
func (curve *KoblitzCurve) ScalarBaseMult(k []byte) (*big.Int, *big.Int) {
	//对输入的大端整数 k 进行模数规约，以确保其在进行椭圆曲线标量乘运算时的合法性和正确性。
	newK := curve.moduloReduce(k)
	//为了确定进行椭圆曲线标量乘运算时需要对输入整数 k 进行的补齐操作数量
	diff := len(curve.bytePoints) - len(newK)

	// Point Q = ∞ (point at infinity).
	qx, qy, qz := new(FieldVal), new(FieldVal), new(FieldVal)

	// curve.bytePoints has all 256 byte points for each 8-bit window. The
	// strategy is to add up the byte points. This is best understood by
	// expressing k in base-256 which it already sort of is.
	// Each "digit" in the 8-bit window can be looked up using bytePoints
	// and added together.
	//通过 newK 中的每个字节 byteVal 来确定基点数组 curve.bytePoints 中的索引。
	/*根据索引从 curve.bytePoints 中获取相应的基点坐标。
	将获取的基点坐标与当前的累加结果 Q 进行 Jacobian 坐标系下的点加运算，结果仍然保存在 Q 中
	循环结束后，将累加的结果 Q 转换为仿射坐标，并作为函数的返回值。*/
	for i, byteVal := range newK {
		//curve.bytePoints 是一个二维切片，存储了预先计算好的椭圆曲线上的点。在 Koblitz 椭圆曲线上，通常会预先计算一些点以加速标量乘运算
		p := curve.bytePoints[diff+i][byteVal]
		curve.addJacobian(qx, qy, qz, &p[0], &p[1], &p[2], qx, qy, qz)
	}
	return curve.fieldJacobianToBigAffine(qx, qy, qz)
}

func (curve *KoblitzCurve) ScalarBaseMultPoint(k []byte) *Point {
	newK := curve.moduloReduce(k)
	diff := len(curve.bytePoints) - len(newK)

	// Point Q = ∞ (point at infinity).
	qx, qy, qz := new(FieldVal), new(FieldVal), new(FieldVal)

	// curve.bytePoints has all 256 byte points for each 8-bit window. The
	// strategy is to add up the byte points. This is best understood by
	// expressing k in base-256 which it already sort of is.
	// Each "digit" in the 8-bit window can be looked up using bytePoints
	// and added together.
	for i, byteVal := range newK {
		p := curve.bytePoints[diff+i][byteVal]
		curve.addJacobian(qx, qy, qz, &p[0], &p[1], &p[2], qx, qy, qz)
	}
	px, py := curve.fieldJacobianToBigAffine(qx, qy, qz)
	return &Point{px, py}
}

func (curve *KoblitzCurve) ScalarBaseMultField(k []byte) (*FieldVal, *FieldVal, *FieldVal) {
	newK := curve.moduloReduce(k)
	diff := len(curve.bytePoints) - len(newK)

	// Point Q = ∞ (point at infinity).
	qx, qy, qz := new(FieldVal), new(FieldVal), new(FieldVal)

	// curve.bytePoints has all 256 byte points for each 8-bit window. The
	// strategy is to add up the byte points. This is best understood by
	// expressing k in base-256 which it already sort of is.
	// Each "digit" in the 8-bit window can be looked up using bytePoints
	// and added together.
	for i, byteVal := range newK {
		p := curve.bytePoints[diff+i][byteVal]
		curve.addJacobian(qx, qy, qz, &p[0], &p[1], &p[2], qx, qy, qz)
	}
	return qx, qy, qz
}

// QPlus1Div4 returns the (P+1)/4 constant for the curve for use in calculating
// square roots via exponentiation.
//
// DEPRECATED: The actual value returned is (P+1)/4, where as the original
// method name implies that this value is (((P+1)/4)+1)/4. This method is kept
// to maintain backwards compatibility of the API. Use Q() instead.
func (curve *KoblitzCurve) QPlus1Div4() *big.Int {
	return curve.q
}

// Q returns the (P+1)/4 constant for the curve for use in calculating square
// roots via exponentiation.
func (curve *KoblitzCurve) Q() *big.Int {
	return curve.q
}

var initonce sync.Once
var secp256k1 KoblitzCurve

func initAll() {
	initS256()
}

// fromHex converts the passed hex string into a big integer pointer and will
// panic is there is an error.  This is only provided for the hard-coded
// constants so errors in the source code can bet detected. It will only (and
// must only) be called for initialization purposes.
func fromHex(s string) *big.Int {
	r, ok := new(big.Int).SetString(s, 16)
	if !ok {
		panic("invalid hex in source file: " + s)
	}
	return r
}

func initS256() {
	// Curve parameters taken from [SECG] section 2.4.1.
	//参数定义了 secp256k1 椭圆曲线 初始化
	secp256k1.CurveParams = new(elliptic.CurveParams)
	secp256k1.P = fromHex("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F")
	secp256k1.N = fromHex("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141")
	secp256k1.B = fromHex("0000000000000000000000000000000000000000000000000000000000000007")
	secp256k1.Gx = fromHex("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798")
	secp256k1.Gy = fromHex("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8")
	secp256k1.BitSize = 256
	secp256k1.q = new(big.Int).Div(new(big.Int).Add(secp256k1.P,
		big.NewInt(1)), big.NewInt(4))
	secp256k1.H = 1
	secp256k1.halfOrder = new(big.Int).Rsh(secp256k1.N, 1)
	secp256k1.fieldB = new(FieldVal).SetByteSlice(secp256k1.B.Bytes())

	// Provided for convenience since this gets computed repeatedly.
	//曲线的字节长度
	secp256k1.byteSize = secp256k1.BitSize / 8

	// Deserialize and set the pre-computed table used to accelerate scalar
	// base multiplication.  This is hard-coded data, so any errors are
	// panics because it means something is wrong in the source code.
	if err := loadS256BytePoints(); err != nil {
		panic(err)
	}

	// Next 6 constants are from Hal Finney's bitcointalk.org post:
	// https://bitcointalk.org/index.php?topic=3238.msg45565#msg45565
	// May he rest in peace.
	//
	// They have also been independently derived from the code in the
	// EndomorphismVectors function in gensecp256k1.go.
	secp256k1.lambda = fromHex("5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72")
	secp256k1.beta = new(FieldVal).SetHex("7AE96A2B657C07106E64479EAC3434E99CF0497512F58995C1396C28719501EE")
	secp256k1.a1 = fromHex("3086D221A7D46BCDE86C90E49284EB15")
	secp256k1.b1 = fromHex("-E4437ED6010E88286F547FA90ABFE4C3")
	secp256k1.a2 = fromHex("114CA50F7A8E2F3F657C1108D9D44CFD8")
	secp256k1.b2 = fromHex("3086D221A7D46BCDE86C90E49284EB15")

	// Alternatively, we can use the parameters below, however, they seem
	//  to be about 8% slower.
	// secp256k1.lambda = fromHex("AC9C52B33FA3CF1F5AD9E3FD77ED9BA4A880B9FC8EC739C2E0CFC810B51283CE")
	// secp256k1.beta = new(FieldVal).SetHex("851695D49A83F8EF919BB86153CBCB16630FB68AED0A766A3EC693D68E6AFA40")
	// secp256k1.a1 = fromHex("E4437ED6010E88286F547FA90ABFE4C3")
	// secp256k1.b1 = fromHex("-3086D221A7D46BCDE86C90E49284EB15")
	// secp256k1.a2 = fromHex("3086D221A7D46BCDE86C90E49284EB15")
	// secp256k1.b2 = fromHex("114CA50F7A8E2F3F657C1108D9D44CFD8")
}

// S256 returns a Curve which implements secp256k1.
func S256() *KoblitzCurve {
	initonce.Do(initAll)
	return &secp256k1
}
