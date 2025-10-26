package btcec

import (
	"bufio"
	"encoding/binary"
	"encoding/gob"
	"fmt"
	"io"
	"math/big"
	"os"
	"runtime"
	"strings"
	"sync"
	"time"
)

const (
	Ilen      = 13    // l2-1 10 多少层 log2N 18
	Imax      = 8192  // 1<< Ilen //1024  262144
	Treelen   = 16384 // imax*2 //2048  524288
	TestNum   = 8192
	Jlen      = 24       // l1-1
	Jmax      = 16777216 // 1<<Jlen
	cuckoolen = 21810381 //jmax*1.3

)

var Threadnum = 4

// var Routinenum = 8
var WG sync.WaitGroup
var IsParTree = 0
var IsNormal = 0

var Mmax = Imax * Jmax * 2
var mflag = int64(Mmax)
var c = S256()

var T1 = Op_NewCuckoo()

var T2x [Imax]FieldVal
var T2y [Imax]FieldVal

var MapT1 map[uint64]uint32

var zero = big.NewInt(0)
var one = big.NewInt(1)
var three = big.NewInt(3)
var seven = big.NewInt(7)

type Point struct {
	X *big.Int
	Y *big.Int
}

type Cipher struct {
	C1x *big.Int
	C1y *big.Int
	c2x *big.Int
	c2y *big.Int
}

type FieldCipher struct {
	c1x *FieldVal
	c1y *FieldVal
	c1z *FieldVal
	c2x *FieldVal
	c2y *FieldVal
	c2z *FieldVal
}

func BuildTree(zs []*FieldVal) [Treelen]*FieldVal {
	var ZTree [Treelen]*FieldVal
	for i := 0; i < Imax; i++ {
		ZTree[i] = zs[i] //初始化叶子结点值
	}
	offset := Imax
	treelen := Imax*2 - 3 //计算树的节点总数
	//treelen1 := treelen - 1
	for i := 0; i < treelen; i += 2 { //构建分支节点
		z := new(FieldVal)
		zmult := z.Mul2(ZTree[i], ZTree[i+1]) //两个节点的父亲节点为两个节点的乘积
		//zmult.Normalize()
		ZTree[offset] = zmult
		offset = offset + 1
	}
	return ZTree
}

func CBuildTree(zs []*FieldVal) []*FieldVal {
	n := Imax

	// 计算叶子节点数，如果是单数则添加一个乘1的节点
	if n%2 != 0 {
		one := new(FieldVal)
		one.SetInt(1) // 假设有一个方法可以将FieldVal设置为1
		zs = append(zs, one)
		n++
	}

	// 初始化树的节点列表，树的节点总数未知，用切片动态构建
	ZTree := make([]*FieldVal, 0)
	ZTree = append(ZTree, zs...)

	currentLevelSize := n

	// 构建树的层次
	for currentLevelSize > 1 {

		if currentLevelSize%2 != 0 {
			one := new(FieldVal)
			one.SetInt(1)
			ZTree = append(ZTree, one)
			currentLevelSize++
		}
		nextLevelStart := len(ZTree) // 下一层的起始位置

		// 计算当前层的父节点并添加到树中
		for i := nextLevelStart - currentLevelSize; i < nextLevelStart; i += 2 {
			z := new(FieldVal)
			zmult := z.Mul2(ZTree[i], ZTree[i+1]) // 两个节点的父亲节点为两个节点的乘积
			ZTree = append(ZTree, zmult)
		}
		currentLevelSize /= 2
	}

	return ZTree
}

func ABuildTree(zs []*FieldVal) ([]*FieldVal, []int) {
	n := len(zs)

	// 计算叶子节点数，如果是单数则添加一个乘1的节点
	if n%2 != 0 {
		one := new(FieldVal)
		one.SetInt(1) // 假设有一个方法可以将FieldVal设置为1
		zs = append(zs, one)
		n++
	}

	// 初始化树的节点列表，树的节点总数未知，用切片动态构建
	ZTree := make([]*FieldVal, 0)
	ZTree = append(ZTree, zs...)

	currentLevelSize := n
	levelSizes := []int{currentLevelSize} // 记录每层节点数的数组

	// 构建树的层次
	for currentLevelSize > 1 {
		if currentLevelSize%2 != 0 {
			one := new(FieldVal)
			one.SetInt(1)
			ZTree = append(ZTree, one)
			currentLevelSize++
		}
		nextLevelStart := len(ZTree) // 下一层的起始位置

		// 计算当前层的父节点并添加到树中
		for i := nextLevelStart - currentLevelSize; i < nextLevelStart; i += 2 {
			z := new(FieldVal)
			zmult := z.Mul2(ZTree[i], ZTree[i+1]) // 两个节点的父亲节点为两个节点的乘积
			ZTree = append(ZTree, zmult)
		}
		currentLevelSize /= 2
		levelSizes = append(levelSizes, currentLevelSize) // 记录当前层节点数
	}

	return ZTree, levelSizes
}

// 修改的big.Int类型的构建树的方法
func kBuildTree(zs []*big.Int) [Treelen]*big.Int {
	var ZTree [Treelen]*big.Int
	for i := 0; i < Imax; i++ {
		ZTree[i] = zs[i] //初始化叶子结点值
	}
	offset := Imax
	treelen := Imax*2 - 3 //计算树的节点总数
	//treelen1 := treelen - 1
	for i := 0; i < treelen; i += 2 { //构建分支节点
		z := new(big.Int)
		zmult := z.Mul(ZTree[i], ZTree[i+1]) //两个节点的父亲节点为两个节点的乘积
		//zmult.Normalize()
		ZTree[offset] = zmult
		offset = offset + 1
	}
	return ZTree
}

func GetInvTree(rootinv *FieldVal, ZTree [Treelen]*FieldVal) [Treelen]*FieldVal { //逆树是从顶向下构建
	var ZinvTree [Treelen]*FieldVal
	treelen := Imax*2 - 2 //逆树的总节点数
	prevfloorflag := treelen
	prevfloornum := 1
	thisfloorflag := treelen
	treeroot_inv := new(FieldVal)
	treeroot_inv.Set(rootinv)
	ZinvTree[prevfloorflag] = treeroot_inv //将根节点的逆元存入
	for i := 0; i < Ilen; i++ {
		thisfloornum := prevfloornum * 2             //计算当前层的节点数，是上一层节点数的两倍。
		thisfloorflag = prevfloorflag - thisfloornum //计算当前层的起始位置
		for f := 0; f < thisfloornum; f++ {
			thisindex := f + thisfloorflag
			ztreeindex := thisindex ^ 1
			z := new(FieldVal)
			ZinvTree[thisindex] = z.Mul2(ZTree[ztreeindex], ZinvTree[prevfloorflag+(f/2)])
		}
		prevfloorflag = thisfloorflag
		prevfloornum = prevfloornum * 2
	}
	return ZinvTree
}

func AGetInvTree(rootinv *FieldVal, ZTree []*FieldVal, levelSizes []int) []*FieldVal {
	treelen := len(ZTree) - 1 // 逆树的总节点数
	ZinvTree := make([]*FieldVal, len(ZTree))

	// 初始化根节点的逆元
	treeroot_inv := new(FieldVal)
	treeroot_inv.Set(rootinv)
	ZinvTree[treelen] = treeroot_inv // 将根节点的逆元存入

	prevfloorflag := treelen // 上一层的起始位置
	for level := len(levelSizes) - 2; level >= 0; level-- {
		thisfloornum := levelSizes[level]             // 当前层的节点数
		thisfloorflag := prevfloorflag - thisfloornum // 当前层的起始位置

		// 处理当前层的节点
		for f := 0; f < thisfloornum; f++ {
			thisindex := f + thisfloorflag
			z := new(FieldVal)
			ZinvTree[thisindex] = z.Mul2(ZinvTree[prevfloorflag+(f/2)], ZTree[thisindex^1])
		}

		prevfloorflag = thisfloorflag
	}

	return ZinvTree
}

func CGetInvTree(rootinv *FieldVal, ZTree []*FieldVal) []*FieldVal { //逆树是从顶向下构建
	treelen := len(ZTree) - 1 // 逆树的总节点数
	ZinvTree := make([]*FieldVal, len(ZTree))

	prevfloorflag := treelen
	prevfloornum := 1
	thisfloorflag := treelen
	treeroot_inv := new(FieldVal)
	treeroot_inv.Set(rootinv)
	ZinvTree[prevfloorflag] = treeroot_inv //将根节点的逆元存入
	// 计算树的层数 Ilen
	high := 0
	temp := len(ZTree)
	for temp > 1 {
		high++
		temp = (temp + 1) / 2
	}
	for i := 0; i < high; i++ {
		thisfloornum := prevfloornum * 2             //计算当前层的节点数，是上一层节点数的两倍。
		thisfloorflag = prevfloorflag - thisfloornum //计算当前层的起始位置
		for f := 0; f < thisfloornum; f++ {
			thisindex := f + thisfloorflag
			ztreeindex := thisindex ^ 1
			z := new(FieldVal)
			ZinvTree[thisindex] = z.Mul2(ZTree[ztreeindex], ZinvTree[prevfloorflag+(f/2)])
		}
		prevfloorflag = thisfloorflag
		prevfloornum = prevfloornum * 2
	}
	return ZinvTree
}

// 修改参数类型 仿射坐标下的构建逆树方法
func kGetInvTree(rootinv *big.Int, ZTree [Treelen]*big.Int) [Treelen]*big.Int { //逆树是从顶向下构建
	var ZinvTree [Treelen]*big.Int
	treelen := Imax*2 - 2 //逆树的总节点数
	prevfloorflag := treelen
	prevfloornum := 1
	thisfloorflag := treelen
	treeroot_inv := new(big.Int)
	treeroot_inv.Set(rootinv)
	ZinvTree[prevfloorflag] = treeroot_inv //将根节点的逆元存入
	for i := 0; i < Ilen; i++ {
		thisfloornum := prevfloornum * 2             //计算当前层的节点数，是上一层节点数的两倍。
		thisfloorflag = prevfloorflag - thisfloornum //计算当前层的起始位置
		for f := 0; f < thisfloornum; f++ {
			thisindex := f + thisfloorflag
			ztreeindex := thisindex ^ 1 //得到二叉树树做对称的位置 好去找正常树要相乘的点
			z := new(big.Int)
			ZinvTree[thisindex] = z.Mul(ZTree[ztreeindex], ZinvTree[prevfloorflag+(f/2)])
		}
		prevfloorflag = thisfloorflag
		prevfloornum = prevfloornum * 2
	}
	return ZinvTree
}

func Encrypt(pubkey *PublicKey, m *big.Int) *Cipher { // 加密函数 传参是公钥 和明文，返回值是密文 c1 和 c2
	start1 := time.Now()
	r, _ := NewPrivateKey(c)                                    //随机生成这个曲线的 加密时的随机数r
	rpkx, rpky := c.ScalarMult(pubkey.X, pubkey.Y, r.D.Bytes()) //rQ Q=kP k是私钥 r是随机数
	cost1 := time.Since(start1)
	//fmt.Printf("btcecc encrypt cost=[%s]\n", cost1)
	GetEnc = GetEnc + cost1                 //记录时间变化
	mGx, mGy := c.ScalarBaseMult(m.Bytes()) //记录mG  secp256k1 曲线的G是预定义公开的
	//fmt.Println(mGx)
	if m.Cmp(zero) == -1 {
		mGy = mGy.Sub(c.P, mGy)
	}
	c2x, c2y := c.Add(mGx, mGy, rpkx, rpky) // c2 是mg + rq
	return &Cipher{r.PublicKey.X, r.PublicKey.Y, c2x, c2y}
}

func NormalEnc(pubkey *PublicKey, m *big.Int) *Cipher {
	start1 := time.Now()
	r, _ := NewPrivateKey(c)
	rpkx, rpky := c.ScalarMult(pubkey.X, pubkey.Y, r.D.Bytes())
	cost1 := time.Since(start1)
	//fmt.Printf("btcecc encrypt cost=[%s]\n", cost1)
	GetEnc = GetEnc + cost1
	mGx, mGy := c.ScalarBaseMult(m.Bytes())
	c2x, c2y := c.Add(mGx, mGy, rpkx, rpky)
	return &Cipher{r.PublicKey.X, r.PublicKey.Y, c2x, c2y}
}

func EncryptJob(pubkey *PublicKey, m *big.Int) *FieldCipher {
	r, _ := NewPrivateKey(c)
	rpkx, rpky := c.ScalarMult(pubkey.X, pubkey.Y, r.D.Bytes())
	mGx, mGy := c.ScalarBaseMult(m.Bytes())
	if m.Cmp(zero) == -1 {
		mGy = mGy.Sub(c.P, mGy)
	}
	c2x, c2y, c2z := c.Add1(mGx, mGy, rpkx, rpky)
	c1x, c1y := c.bigAffineToField(r.PublicKey.X, r.PublicKey.Y)
	c1z := new(FieldVal).SetInt(1)
	return &FieldCipher{c1x, c1y, c1z, c2x, c2y, c2z}
}

var GetEnc time.Duration = 0
var GetmG time.Duration = 0
var GetX21 time.Duration = 0
var GetTree1 time.Duration = 0
var GetInv time.Duration = 0
var GetTree2 time.Duration = 0
var BSGS time.Duration = 0
var Verify time.Duration = 0
var GetHash time.Duration = 0
var GetX3 time.Duration = 0
var GetSearch time.Duration = 0

func BytesToUint32(b []byte) Hash {
	_ = b[3] //这是一个占位符语句，它确保字节数组 b 的长度至少为4
	return Hash(b[3]) | Hash(b[2])<<8 | Hash(b[1])<<16 | Hash(b[0])<<24
}

func (curve *KoblitzCurve) kGetx3y3(x1, y1, x2, y2, x3, y3 *FieldVal) {
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

// 批量获得两不同点相加
func Get3(x1, y1, x2, y2, x3, y3 []*big.Int, ZinvTree [Treelen]*big.Int, start, end int, overthreadnum chan int) {
	//time.Sleep(time.Duration(4) * time.Second)

	for j := start; j < end; j++ {
		var h, i, r, v, d, f, g big.Int
		if x3[j] != nil && y3[j] != nil {
			continue
		}
		r.Sub(y2[j], y1[j])          //r = y2 - y1
		h.Mul(&r, ZinvTree[j])       //h = (y2 - y1)/(x2 - x1)
		h.Mod(&h, c.P)               //h = (y2 - y1)/(x2 - x1) mod p
		i.Mul(&h, &h)                // h = h^2
		v.Add(x1[j], x2[j]).Neg(&v)  //v = -x1 - x2
		f.Set(&i).Add(&i, &v)        // f = h^2 - x1 - x2
		x3[j].Mod(&f, c.P)           // x3 = h^2 - x1 - x2 mod p
		d.Sub(x1[j], x3[j])          //d = x1 -x3
		g.Mul(&d, &h).Sub(&g, y1[j]) // g = h*(x1 -x3) -y1
		y3[j].Mod(&g, c.P)           //y3 = h*(x1 -x3) -y1 mod p

	}
	//overthreadnum <- 1
}

// 批量获得两不同点相加
func fGet3(x1, y1, x2, y2, x3, y3 []*FieldVal, ZinvTree [Treelen]*FieldVal, start, end int, overthreadnum chan int) {
	//time.Sleep(time.Duration(4) * time.Second)

	for j := start; j < end; j++ {
		var h, i, r, v, d, e, l, m FieldVal
		if ZinvTree[j].ToBigInt().Cmp(big.NewInt(1)) == 0 {
			continue
		}

		//fmt.Println("-------------------------------")
		//fmt.Println("ZinvTree[j], x1[j], y1[j], x2[j], y2[j]")
		// g.Set(x1[j]).Negate(1)
		// f.Set(&g).Add(x2[j])
		// fmt.Println("f.Inverse()", f.Inverse())
		//ZinvTree[j]有问题
		//fmt.Println(ZinvTree[j], x1[j], y1[j], x2[j], y2[j])
		r.Set(y1[j]).Negate(1).Add(y2[j]) //r = y2 - y1
		//fmt.Println("r.Set(y1[j].Negate(1)).Add(y2[j])", r)
		h.Mul2(&r, ZinvTree[j]) //h = (y2 - y1)/(x2 - x1)
		//fmt.Println("h.Mul2(&r, ZinvTree[j]) ", h)
		h.Normalize()
		//fmt.Println("h.Normalize() ", h)
		i.SquareVal(&h) // h = h^2
		//fmt.Println("i.SquareVal(&h) ", i)
		//i.Normalize()
		//fmt.Println("i.Normalize() ", i)
		e.Set(x1[j]).Add(x2[j])
		//fmt.Println("e.Set(x1[j]).Add(x2[j])", e)
		//e.Normalize()
		//fmt.Println("e.Normalize()", e)
		v.Set(&e).Negate(1) //v = -x1 - x2
		//fmt.Println("v.Set(e.Negate(1)) ", v)
		l.Set(&i).Add(&v) // x3 = h^2 - x1 - x2
		//fmt.Println("l.Set(&i).Add(&v) ", l)
		l.Normalize()
		//fmt.Println("l.Normalize() ", l)

		d.Set(&l).Negate(1).Add(x1[j]) //d = x1 -x3
		//fmt.Println("d.Set(&l).Negate(1).Add(x1[j]) ", d)
		m.Set(y1[j]).Negate(1)
		//fmt.Println("m.Set(y1[j]).Negate(1) ", m)
		//m.Normalize()
		//fmt.Println("m.Normalize() ", m)
		y3[j].Mul2(&h, &d).Add(&m) // g= h*(x1 -x3) -y1
		//fmt.Println("y3[j].Mul2(&h, &d).Add((y1[j]).Negate(1)) ", y3[j])
		y3[j].Normalize()
		//fmt.Println("y3[j].Normalize() ", y3[j])
		x3[j].Set(&l)
		//fmt.Println("x3[j].Set(&l) ", x3[j])
		//x3[j].Normalize()
		//fmt.Println("x3[j].Normalize() ", x3[j])

	}
	//overthreadnum <- 1
}

// spec256k1 a = 0 b = 7  y^2 = x^3 + 7
// 批量获得不同点的倍乘
func Get2(x1, y1, x2, y2 []*big.Int, ZinvTree [Treelen]*big.Int, start, end int, overthreadnum chan int) {
	//time.Sleep(time.Duration(4) * time.Second)

	if x1[0] == nil { //假如x1为空指针就执行一下操作
		x2 = x1
		y2 = y1
		return
	}

	for j := start; j < end; j++ {
		var h, i, r, v, d, f, g big.Int
		if x2[j] != nil && y2[j] != nil {
			continue
		}

		r.Mul(x1[j], x1[j]).Mul(&r, big.NewInt(3)) // r = x1^2 * 3 a = 0

		h.Mul(&r, ZinvTree[j])       // h = (x1^2 * 3 + a)/ 2y1
		h.Mod(&h, c.P)               // h = (x1^2 * 3 + a)/ 2y1 mod p
		i.Mul(&h, &h)                // h = h^2
		v.Add(x1[j], x1[j]).Neg(&v)  //v = -x1 - x1
		f.Set(&i).Add(&i, &v)        // f = h^2 - x1 - x1
		x2[j].Mod(&f, c.P)           // x2 = h^2 - x1 - x1 mod p
		d.Sub(x1[j], x2[j])          //d = x1 -x2
		g.Mul(&d, &h).Sub(&g, y1[j]) // g = h*(x1 -x2) -y1
		y2[j].Mod(&g, c.P)           //y2 = h*(x1 -x2) -y1 mod p

	}

	//overthreadnum <- 1
}

// spec256k1 a = 0 b = 7  y^2 = x^3 + 7
// 批量获得不同点的倍乘 FieldVal类型
func fGet2(x1, y1, x2, y2 []*FieldVal, ZinvTree [Treelen]*FieldVal, start, end int, overthreadnum chan int) {
	//time.Sleep(time.Duration(4) * time.Second)

	for j := start; j < end; j++ {
		var h, i, r, v, d, k, l, m FieldVal
		// if !x2[j].IsZero() && !y2[j].IsZero() {
		// 	continue
		// }
		if ZinvTree[j].ToBigInt().Cmp(big.NewInt(1)) == 0 {
			continue
		}
		// fmt.Println("-------------------------------")
		// fmt.Println("ZinvTree[j],x1[j],y1[j]")
		// fmt.Println(ZinvTree[j], x1[j], y1[j])
		r.SquareVal(x1[j]).MulInt(3) // r = x1^2 * 3 a = 0
		// fmt.Println("r.SquareVal(x1[j]).MulInt(3)", r)
		h.Mul2(&r, ZinvTree[j]) // h = (x1^2 * 3 + a)/ 2y1
		//h.Set(&r).Mul(ZinvTree[j])
		//fmt.Println("h.Mul2(&r, ZinvTree[j])", h)
		//h.Normalize()
		//fmt.Println("h.Normalize()", h)
		i.SquareVal(&h) // h = h^2
		//fmt.Println("i.SquareVal(&h)", i)
		//i.Normalize()
		//fmt.Println("i.Normalize()", i)
		//************这里x1[j]被改变了  巨坑**********o.Set(x1[j].MulInt(2))
		//o.Set(x1[j].MulInt(2))
		//o.Set(x1[j]).MulInt(2)
		//fmt.Println("o.Set(x1[j].MulInt(2))", o)
		//o.Normalize()
		//fmt.Println("o.Normalize()", o)
		//v.Set(x1[j]).MulInt(2).Negate(1) //v = -x1 - x1
		v.Set(x1[j]).Negate(1)
		//fmt.Println("v.Set(o.Negate(1))", v)
		//v.Normalize()
		//fmt.Println("v.Normalize()", v)
		//fmt.Println("x1[j]", x1[j])
		//这里把x1的值改变了 x2[j].Set(&i).Add(&v) // x2 = h^2 - x1 - x1
		//l.Set(&i).Add(&v) // x2 = h^2 - x1 - x1
		l.Set(&i).Add(&v)
		l.Add(&v)
		//fmt.Println("x1[j]", x1[j])
		//fmt.Println("l.Set(&i).Add(&v)", l)
		l.Normalize()
		//fmt.Println("l.Normalize()", l)
		//这里有问题 涉及到负数就要及时Normalize
		k.Set(&l).Negate(1) //-x2
		//fmt.Println("k.Set(x2[j].Negate(1))", k)
		//k.Normalize()
		//fmt.Println("k.Normalize()", k)
		// d.Set(&k).Add(x1[j]) //d = x1 -x2
		//fmt.Println("x1[j]", x1[j])
		d.Set(&k).Add(x1[j]) //d = x1 -x2
		//fmt.Println("d.Set(&k).Add(x1[j])", d)
		//d.Normalize()
		//fmt.Println("d.Normalize()", d)
		//g.Mul2(&h, &d) // g = h*(x1 -x2)
		m.Set(y1[j]).Negate(1)
		// m.Normalize()
		y2[j].Mul2(&h, &d) // g = h*(x1 -x2)
		//fmt.Println("g.Mul2(&h, &d)", g)
		//y2[j].Set(y1[j]).Negate(1).Add(&g) // y2 = h*(x1 -x2) -y1
		y2[j].Add(&m)

		// y2[j].Add(&m)
		//fmt.Println("y2[j].Set(y1[j].Negate(1)).Add(&g)", y2[j])
		y2[j].Normalize()
		//fmt.Println("y2[j].Normalize()", y2[j])
		x2[j].Set(&l)
		//fmt.Println("x2[j].Set(&l)", x2[j])
		//x2[j].Normalize()
		//fmt.Println("x2[j].Normalize()", x2[j])

	}

	//overthreadnum <- 1
}

func FGet2(x1, y1, x2, y2 []*FieldVal, ZinvTree [Treelen]*FieldVal, start, end int, overthreadnum chan int) {
	//time.Sleep(time.Duration(4) * time.Second)

	// if x1[0] == nil { //假如x1为空指针就执行一下操作
	// 	x2 = x1
	// 	y2 = y1
	// 	return
	// }

	for j := start; j < end; j++ {
		var h, i, r, v, d, g FieldVal
		// if !x2[j].IsZero() && !y2[j].IsZero() {
		// 	continue
		// }
		if ZinvTree[j].ToBigInt().Cmp(big.NewInt(1)) == 0 {
			continue
		}

		r.SquareVal(x1[j]).MulInt(3) // r = x1^2 * 3 a = 0
		h.Mul2(&r, ZinvTree[j])      // h = (x1^2 * 3 + a)/ 2y1
		h.Normalize()
		i.SquareVal(&h)                  // h = h^2
		v.Set(x1[j]).MulInt(2).Negate(1) //v = -x1 - x1
		x2[j].Set(&i).Add(&v)            // x2 = h^2 - x1 - x1
		x2[j].Normalize()
		d.Set(x2[j].Negate(1)).Add(x1[j])  //d = x1 -x2
		g.Mul2(&h, &d)                     // g = h*(x1 -x2)
		y2[j].Set(y1[j].Negate(1)).Add(&g) // y2 = h*(x1 -x2) -y1
		y2[j].Normalize()

	}

	//overthreadnum <- 1
}

// *********************trick的核心函数FastECDLP***********************
func GetM(mGx *big.Int, ZinvTree [Treelen]*FieldVal, p *FieldVal, fmGx, fmGy *FieldVal, t, start, end, jmax int, m_new []int64, m_bool []bool, overthreadnum chan int) {
	//time.Sleep(time.Duration(4) * time.Second)
	for j := start; j < end; j++ {
		if j == 0 {
			leftxbytes := mGx.Bytes()
			i, ok := T1.Op_search(leftxbytes) //假如mg有和T1相等的 直接的得出m并验证 在表存寻找是否有相同的秘钥
			if ok {
				m := int64(i)
				TestmGx, _ := c.ScalarBaseMult(big.NewInt(m).Bytes())
				r1 := mGx.Cmp(TestmGx)
				if r1 == 0 {
					m_bool[t] = true
					m_new[t] = m
					break
				}
			}
		}

		ft2x, ft2y := T2x[j], T2y[j]
		/*
			ft2x := T2x[j]
			t2x := new(big.Int).SetBytes(ft2x.Bytes()[:])
			t2y := GetY(t2x)
			ft2y := new(FieldVal).SetByteSlice(t2y.Bytes())
		*/

		leftx, invleftx := c.NewGetx3(fmGx, fmGy, &ft2x, &ft2y, ZinvTree[j], p) //求出mG与T2 的加法的结果  得出一个正常结果 和一个求负结果
		leftxbytes := leftx.Bytes()

		i, ok := T1.Op_search(leftxbytes) //然后再T1中找和 上面得出结果是否有一样的 一样的话就按照公式得出结果
		if ok {
			m1 := int64(j+1) * int64(jmax) * 2
			m2 := int64(i)
			m := m1 + m2
			TestmGx, _ := c.ScalarBaseMult(big.NewInt(m).Bytes())
			r1 := mGx.Cmp(TestmGx)
			if r1 == 0 {
				m_bool[t] = true
				m_new[t] = m
				break
			}

			m = m1 - m2
			TestmGx, _ = c.ScalarBaseMult(big.NewInt(m).Bytes())
			r1 = mGx.Cmp(TestmGx)
			if r1 == 0 {
				m_bool[t] = true
				m_new[t] = m
				break
			}
		}
		leftxbytes = invleftx.Bytes()
		i, ok = T1.Op_search(leftxbytes)
		if ok {
			m1 := int64(j+1) * int64(jmax) * 2
			m2 := int64(i)
			m := -(m1 + m2)
			TestmGx, _ := c.ScalarBaseMult(big.NewInt(m).Bytes())
			r1 := mGx.Cmp(TestmGx)
			if r1 == 0 {
				m_bool[t] = true
				m_new[t] = m
				break
			}

			m = -(m1 - m2)
			TestmGx, _ = c.ScalarBaseMult(big.NewInt(m).Bytes())
			r1 = mGx.Cmp(TestmGx)
			if r1 == 0 {
				m_bool[t] = true
				m_new[t] = m
				break
			}
		}
	}
	overthreadnum <- 1
}

func ParDecrypt(priv *PrivateKey, cipher *Cipher) (*big.Int, string) { //传入参数是秘钥 和密文 c1 c2

	//start1 := time.Now()
	var m int64 = mflag
	skc1x, skc1y := c.ScalarMult(cipher.C1x, cipher.C1y, priv.D.Bytes()) //标量乘的做法， 求私钥。 D为随机的大整数， c1为随机数乘G点， skc1 算的是sk*c1的值
	if skc1x.Cmp(cipher.c2x) == 0 {                                      //比较 因为m*G = c2 - skc1. 假如他们相当的话 明文就是zero
		return zero, ""
	}
	//ZTree := make([]*FieldVal, Imax*2, Imax*2)
	//ZinvTree := make([]*FieldVal, Imax*2, Imax*2)
	//var ZTree [Treelen]*FieldVal
	//var ZinvTree [Treelen]*FieldVal
	inv_skc1y := new(big.Int)
	inv_skc1y.Add(c.P, inv_skc1y)                               //p是循环群p
	inv_skc1y.Sub(inv_skc1y, skc1y)                             //这几步是算-y
	mGx, mGy := c.Add(cipher.c2x, cipher.c2y, skc1x, inv_skc1y) //然后计算c2 - sk*c1 实际上计算c2 + （x，-y）
	//fmt.Println(mGx)
	fmGx, fmGy := c.bigAffineToField(mGx, mGy) //把坐标系上的点转换到有限域上
	zs := make([]*FieldVal, Imax)              //创建一个切片 长度为Imax
	//直接寻找在T2上和mG是否有相同的点 有相同的话 直接得出明文
	for i := 0; i < Imax; i++ {
		ft2x := T2x[i]               //把T2的x值赋值
		zs[i] = c.Getz3(fmGx, &ft2x) //计算H ft2x - fmGx
		zs[i].Normalize()
		if zs[i].Equals(fieldZero) == true { //如果zs[i]的值为零 说明有i的时候mG和T2相等 就直接解密得出m
			m = int64(Jmax*2) * int64(i+1)
			mbigint := big.NewInt(m)
			_, TestmGy := c.ScalarBaseMult(mbigint.Bytes()) //标量乘m*G
			r1 := mGx.Cmp(TestmGy)                          //测试mg和算出来的m *g  测试m是否正确
			if r1 == 0 {
				WG.Done() //并发程序最后的同步
				return big.NewInt(m), "secuess"
			}
			WG.Done()
			return new(big.Int).Neg(mbigint), "secuess"
		}
	}

	//runtime.GOMAXPROCS(runtime.NumCPU()) trick构建树的操作
	runtime.GOMAXPROCS(Threadnum)
	ZTree := BuildTree(zs)                                        //构建树 是利用大步的T2中与GM x的差来构建 因为求第三点的时候需要用到 x之差的逆
	treeroot_inv := new(FieldVal).Set(ZTree[Treelen-2]).Inverse() //获取树的根节点，然后计算其逆元
	ZinvTree := GetInvTree(treeroot_inv, ZTree)                   //构建逆树
	p := new(FieldVal).SetByteSlice(c.P.Bytes())                  //将p点构建转化为有限域上

	overthreadnum := make(chan int, Threadnum) //创建和线程数相同的缓冲通道
	batch := Imax / (Threadnum)                //分组数
	m_new := make([]int64, Threadnum)
	m_bool := make([]bool, Threadnum)
	for t := 0; t < Threadnum; t++ {
		m_new[t] = m
		m_bool[t] = false
	}

	for t := 0; t < Threadnum; t++ { //t循环迭代的编号， batch循环迭代的范围
		go GetM(mGx, ZinvTree, p, fmGx, fmGy, t, t*batch, (t+1)*batch, Jmax, m_new, m_bool, overthreadnum) //得到m  这里用的是BSGS算法 并且用了T1 T2加速计算的 trick
	}

	for i := 0; i < Threadnum; i++ {
		<-overthreadnum
	}
	acc_c := 0
	for t := 0; t < Threadnum; t++ { //检测计算的正确性
		if m_bool[t] {
			acc_c += 1
			m = m_new[t]
			TestmGx, _ := c.ScalarBaseMult(big.NewInt(m).Bytes())
			r1 := mGx.Cmp(TestmGx)
			if r1 == 0 {
				if acc_c > 1 {
					fmt.Println("getM 多次", acc_c)
				}
				WG.Done()
				return big.NewInt(m), "sucess"
			}

		}
	}
	WG.Done()
	fmt.Println("解密失败", acc_c) //直接用acc_c来计数解密失败的次数 因为成功的时候直接返回了 到了这一步的话 acc_c记录的都是解密失败的次数
	return big.NewInt(0), "decrypt error 2"
}

func GetY(t2x *big.Int) *big.Int {
	t2y2 := new(big.Int).Exp(t2x, three, c.P)
	t2y2 = t2y2.Add(t2y2, seven)
	t2y2 = t2y2.Mod(t2y2, c.P)
	t2y := t2y2.Sqrt(t2y2)
	t2y = t2y.Mod(t2y, c.P)
	inv_t2y := t2y.Sub(c.P, t2y)
	return inv_t2y
}

func NormalDecrypt(priv *PrivateKey, cipher *Cipher) (*big.Int, string) {
	var m int64 = mflag
	skc1x, skc1y := c.ScalarMult(cipher.C1x, cipher.C1y, priv.D.Bytes())
	if skc1x.Cmp(cipher.c2x) == 0 {
		return zero, ""
	}
	inv_skc1y := new(big.Int)
	inv_skc1y.Add(c.P, inv_skc1y)
	inv_skc1y.Sub(inv_skc1y, skc1y)
	mGx, mGy := c.Add(cipher.c2x, cipher.c2y, skc1x, inv_skc1y)
	start := time.Now()
	for j := 0; j < Imax; j++ {
		if j == 0 {
			// hash time
			leftxbytes := mGx.Bytes()[:8]
			x64 := binary.BigEndian.Uint64(leftxbytes)
			i, ok := MapT1[x64]
			if ok {
				m = int64(i)
				cost := time.Since(start)
				GetmG = GetmG + cost
				break
			}
		}
		//ft2x, ft2y := T2x[j], T2y[j]
		ft2x := T2x[j]
		t2x := new(big.Int).SetBytes(ft2x.Bytes()[:])
		t2y := GetY(t2x)
		//ft2y := new(big.Int).SetBytes(ft2y.Bytes()[:])
		leftx, _ := c.Add(mGx, mGy, t2x, t2y)
		//leftx := c.Getx3(fmGx, fmGy, ft2x, ft2y, ZinvTree[j])
		leftxbytes := leftx.Bytes()[:8]
		x64 := binary.BigEndian.Uint64(leftxbytes)
		i, ok := MapT1[x64]
		if ok {
			m = int64(j)*int64(Jmax) + int64(i)
			cost := time.Since(start)
			GetmG = GetmG + cost
			break
		}
	}
	return big.NewInt(m), "sucess"
}

func GetZS(c *KoblitzCurve, zs []*FieldVal, fmGx *FieldVal, start int, end int) {
	for i := start; i < end; i++ {
		ft2x := T2x[i]
		zs[i] = c.Getz3(fmGx, &ft2x)
	}
}

func HomoAddField(c1 *Cipher, c2 *Cipher) *FieldCipher {
	c1x, c1y, c1z := c.Add1(c1.C1x, c1.C1y, c2.C1x, c2.C1y)
	c2x, c2y, c2z := c.Add1(c1.c2x, c1.c2y, c2.c2x, c2.c2y)
	return &FieldCipher{c1x, c1y, c1z, c2x, c2y, c2z}
}

func HomoAddField1(c1 *FieldCipher, c2 *FieldCipher) *FieldCipher {
	c1x, c1y, c1z := new(FieldVal), new(FieldVal), new(FieldVal)
	c2x, c2y, c2z := new(FieldVal), new(FieldVal), new(FieldVal)
	c.AddGeneric(c1.c1x, c1.c1y, c1.c1z, c2.c1x, c2.c1y, c2.c1z, c1x, c1y, c1z)
	c.AddGeneric(c1.c2x, c1.c2y, c1.c2z, c2.c2x, c2.c2y, c2.c2z, c2x, c2y, c2z)
	return &FieldCipher{c1x, c1y, c1z, c2x, c2y, c2z}
}

func HomoAdd(c1 *Cipher, c2 *Cipher) *Cipher {
	c1x, c1y := c.Add(c1.C1x, c1.C1y, c2.C1x, c2.C1y)
	c2x, c2y := c.Add(c1.c2x, c1.c2y, c2.c2x, c2.c2y)
	return &Cipher{c1x, c1y, c2x, c2y}
}

func HomoAddPlainText(c1 *Cipher, c2 *big.Int) *Cipher {
	c2x, c2y := c.ScalarBaseMult(c2.Bytes())
	c2x, c2y = c.Add(c1.c2x, c1.c2y, c2x, c2y)
	return &Cipher{c1.C1x, c1.C1y, c2x, c2y}
}

func HomoMul(c1 *Cipher, k *big.Int) *Cipher {
	c1x, c1y := c.ScalarMult(c1.C1x, c1.C1y, k.Bytes())
	c2x, c2y := c.ScalarMult(c1.c2x, c1.c2y, k.Bytes())
	return &Cipher{c1x, c1y, c2x, c2y}
}

func HomoMulField(c1 *Cipher, k *big.Int) *FieldCipher {
	c1x, c1y, c1z := c.ScalarMultField(c1.C1x, c1.C1y, k.Bytes())
	c2x, c2y, c2z := c.ScalarMultField(c1.c2x, c1.c2y, k.Bytes())
	return &FieldCipher{c1x, c1y, c1z, c2x, c2y, c2z}
}

func ConvertCipher(fieldc *FieldCipher) *Cipher {
	c1x, c1y := c.fieldJacobianToBigAffine(fieldc.c1x, fieldc.c1y, fieldc.c1z)
	c2x, c2y := c.fieldJacobianToBigAffine(fieldc.c2x, fieldc.c2y, fieldc.c2z)
	return &Cipher{c1x, c1y, c2x, c2y}
}

func ReadT1() {
	T1_file, _ := os.Open("/home/lgw/go/src/github.com/cuckoobtcec/genlist/Tx24.bin")
	T1_dec := gob.NewDecoder(T1_file)
	T1_dec.Decode(T1)
	T1_file.Close()
}

func ReadT2() {
	var j int64 = 1
	t1lastx, t1lasty := c.ScalarMult(c.Gx, c.Gy, big.NewInt(int64(Jmax)).Bytes())
	t2x, t2y := c.ScalarMult(t1lastx, t1lasty, zero.Bytes())
	for ; j < int64(Imax); j++ {
		if j >= 1 {
			t2x, t2y = c.Add(t2x, t2y, t1lastx, t1lasty)
		}
		inv_t2y := new(big.Int)
		inv_t2y.Add(c.P, inv_t2y)
		inv_t2y.Sub(inv_t2y, t2y)
		ft2x, ft2y := c.bigAffineToField(t2x, inv_t2y)
		T2x[j] = *ft2x
		T2y[j] = *ft2y
	}
}

func ReadT1AsMap() {
	var i int64 = 1
	filename := "/home/lgw/go/src/github.com/cuckoobtcec/genlist/Tx28.txt"
	file, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	x := big.NewInt(0)
	MapT1 = make(map[uint64]uint32)
	rd := bufio.NewReader(file)
	for {
		line, err := rd.ReadString('\n')
		if err != nil || io.EOF == err {
			break
		} else {
			line = strings.Replace(line, "\n", "", -1)
			x, _ = new(big.Int).SetString(line, 10)
			MapT1[x.Uint64()] = uint32(i)
			if i == int64(Jmax) {
				file.Close()
				break
			}
			i++
		}
	}
}

func PathExists(path string) (bool, error) {
	_, err := os.Stat(path)
	if err == nil {
		return true, nil
	}
	if os.IsNotExist(err) {
		return false, nil
	}
	return false, err
}

func ReadT1frombin() {
	path := "../genlist/T1.bin"
	isexist, _ := PathExists(path)
	if isexist == true {
		T1.Load(path)
	}
}

// func init() { //初试函数是读取T1 和 生成枚举T2 并且存储的是T2x  和 -T2y
// 	if IsNormal == 1 { //读取T1
// 		ReadT1AsMap()
// 	} else {
// 		ReadT1frombin()
// 	}
// 	var j int64 = 0
// 	t1lastx, t1lasty := c.ScalarMult(c.Gx, c.Gy, big.NewInt(int64(Jmax*2)).Bytes()) //计算T1的最后点的值 因为T1是小步的枚举 所以直接用G点和最大值标量乘
// 	t2x, t2y := c.ScalarMult(t1lastx, t1lasty, one.Bytes())                         //T2的第一个值和T1末尾是相邻的
// 	for ; j < int64(Imax); j++ {
// 		//fmt.Printf("%d\n", j)
// 		if j >= 1 {
// 			t2x, t2y = c.Add(t2x, t2y, t1lastx, t1lasty)
// 		}
// 		inv_t2y := new(big.Int) //这几步操作是在算-t2y
// 		inv_t2y.Add(c.P, inv_t2y)
// 		inv_t2y.Sub(inv_t2y, t2y)
// 		ft2x, ft2y := c.bigAffineToField(t2x, inv_t2y)
// 		T2x[j] = *ft2x
// 		T2y[j] = *ft2y
// 	}

// }
