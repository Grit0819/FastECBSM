package main

import (
	//"crypto/ecdsa"
	"crypto/rand"
	"crypto/sha256"
	"fmt"
	"log"
	"math/big"
	"runtime"
	"time"

	btcec "github.com/cuckoobtcec"
)

var totalcost3 time.Duration = 0
var totalcost4 time.Duration = 0
var ecctimeall time.Duration = 0
var eccenall time.Duration = 0
var eccdeall time.Duration = 0
var ecchaddall time.Duration = 0
var ecchaddfieldall time.Duration = 0
var ecchmulall time.Duration = 0
var costM time.Duration = 0
var costW time.Duration = 0

var c = btcec.S256()

func testparbtcec(messages [btcec.TestNum]*big.Int) {
	privKey, _ := btcec.NewPrivateKey(c)
	pubKey := privKey.PubKey()
	var ciphers [btcec.TestNum]*btcec.Cipher
	var c1xs, c1ys []*big.Int
	runtime.GOMAXPROCS(btcec.Threadnum)
	for i := 0; i < btcec.TestNum; i++ {
		start1 := time.Now()
		ciphers[i] = btcec.Encrypt(pubKey, messages[i])

		cost1 := time.Since(start1)
		eccenall = eccenall + cost1
	}

	start2 := time.Now()
	for i := 0; i < btcec.TestNum; i++ {
		btcec.WG.Add(1)
		go btcec.ParDecrypt(privKey, ciphers[i])

	}

	btcec.WG.Wait()
	cost2 := time.Since(start2)
	eccdeall = eccenall + cost2

	for i := 0; i < btcec.TestNum; i++ {
		c1xs[i] = ciphers[i].C1x
		c1ys[i] = ciphers[i].C1y
	}

	startM := time.Now()
	c.ScalarMultBatch(c1xs, c1ys, privKey.D.Bytes())

	btcec.WG.Wait()
	cost9 := time.Since(startM)
	costM = cost9

	startW := time.Now()
	for i := 0; i < btcec.TestNum; i++ {
		btcec.WG.Add(1)
		go c.ScalarMult(c1xs[i], c1ys[i], privKey.D.Bytes())
	}

	btcec.WG.Wait()
	cost10 := time.Since(startW)
	costW = cost10

	for i := 0; i < btcec.TestNum; i++ {
		cipher1 := btcec.Encrypt(pubKey, messages[btcec.TestNum-1-i])

		start4 := time.Now()
		btcec.HomoAdd(ciphers[i], cipher1)
		cost4 := time.Since(start4)

		start5 := time.Now()
		btcec.HomoAddField(ciphers[i], cipher1)
		cost5 := time.Since(start5)

		start6 := time.Now()
		btcec.HomoMul(ciphers[i], messages[btcec.TestNum-1-i])
		cost6 := time.Since(start6)

		ecchaddall = ecchaddall + cost4
		ecchaddfieldall = ecchaddfieldall + cost5
		ecchmulall = ecchmulall + cost6
	}
}

func Sign(priv *btcec.PrivateKey, message string) (r, s *big.Int, err error) {
	hash := sha256.Sum256([]byte(message))

	var k, invK *big.Int
	var kx *big.Int
	for {
		k, err = rand.Int(rand.Reader, btcec.S256().N)
		if err != nil {
			return nil, nil, err
		}

		kx, _ = c.ScalarBaseMult(k.Bytes())
		if kx.Sign() == 0 {
			continue
		}

		r = new(big.Int).Mod(kx, btcec.S256().N)
		if r.Sign() == 0 {
			continue
		}

		invK = new(big.Int).ModInverse(k, btcec.S256().N)
		if invK == nil {
			continue
		}

		z := new(big.Int).SetBytes(hash[:])
		s = new(big.Int).Mul(r, priv.D)
		s.Add(s, z)
		s.Mul(s, invK)
		s.Mod(s, btcec.S256().N)
		if s.Sign() != 0 {
			break
		}
	}

	return r, s, nil
}

func Signbatch(priv []*btcec.PrivateKey, message string) (rs, ss []*big.Int, err error) {
	hash := sha256.Sum256([]byte(message))
	ks := make([]*big.Int, btcec.Imax)
	ksBytes := make([][]byte, btcec.Imax)
	rs = make([]*big.Int, btcec.Imax)
	ss = make([]*big.Int, btcec.Imax)
	for i := 0; i < btcec.Imax; i++ {
		rs[i] = new(big.Int)
		ss[i] = new(big.Int)
		ks[i] = new(big.Int)
		var k *big.Int
		k, err = rand.Int(rand.Reader, btcec.S256().N)
		if err != nil {
			return nil, nil, err
		}
		ks[i] = k
		ksBytes[i] = ks[i].Bytes()
	}
	kxs, _ := c.ScalarBaseMultBatchFieldVal(ksBytes)
	for i := 0; i < btcec.Imax; i++ {

		for {
			if kxs[i].Sign() == 0 {
				var k *big.Int
				k, err = rand.Int(rand.Reader, btcec.S256().N)
				if err != nil {
					return nil, nil, err
				}
				kxs[i], _ = c.ScalarBaseMult(k.Bytes())
				continue
			}
			r := new(big.Int).Mod(kxs[i], btcec.S256().N)
			if r.Sign() == 0 {
				continue
			}
			rs[i] = r

			// Calculate s = k^(-1) * (hash + r * privKey.D) mod N
			invK := new(big.Int).ModInverse(ks[i], btcec.S256().N)
			if invK == nil {
				continue
			}

			z := new(big.Int).SetBytes(hash[:])
			s := new(big.Int).Mul(r, priv[i].D)
			s.Add(s, z)
			s.Mul(s, invK)
			s.Mod(s, btcec.S256().N)
			if s.Sign() != 0 {
				ss[i] = s
				break
			}
		}

	}
	return rs, ss, nil
}

func VerifySignatureBatch(Bx, By []*big.Int, message string, rs, ss []*big.Int) bool {

	hash := sha256.Sum256([]byte(message))
	z := new(big.Int).SetBytes(hash[:])
	ws := make([]*big.Int, btcec.Imax)
	u1s := make([]*big.Int, btcec.Imax)
	u2s := make([]*big.Int, btcec.Imax)
	u1sBytes := make([][]byte, btcec.Imax)
	u2sBytes := make([][]byte, btcec.Imax)
	for i := 0; i < btcec.Imax; i++ {

		if rs[i].Sign() <= 0 || rs[i].Cmp(btcec.S256().N) >= 0 {
			return false
		}
		if ss[i].Sign() <= 0 || ss[i].Cmp(btcec.S256().N) >= 0 {
			return false
		}
		ws[i] = new(big.Int)

		ws[i] = new(big.Int).ModInverse(ss[i], btcec.S256().N)
		if ws[i] == nil {
			return false
		}
		u1s[i] = new(big.Int)
		u2s[i] = new(big.Int)

		u1s[i] = new(big.Int).Mul(z, ws[i])
		u1s[i].Mod(u1s[i], btcec.S256().N)
		u2s[i] = new(big.Int).Mul(rs[i], ws[i])
		u2s[i].Mod(u2s[i], btcec.S256().N)
		u1sBytes[i] = u1s[i].Bytes()
		u2sBytes[i] = u2s[i].Bytes()
	}

	x1s, y1s := c.ScalarBaseMultBatchFieldVal(u1sBytes)
	x2s, y2s := c.ScalarMultBatchFieldValBatchQ(Bx, By, u2sBytes)

	for i := 0; i < btcec.Imax; i++ {
		x, y := c.Add(x1s[i], y1s[i], x2s[i], y2s[i])

		if x.Sign() == 0 && y.Sign() == 0 {
			return false
		}

		v := new(big.Int).Mod(x, btcec.S256().N)
		if v.Cmp(rs[i]) != 0 {
			return false
		}

	}

	return true

}

func VerifySignature(pubkey *btcec.PublicKey, message string, r, s *big.Int) bool {
	if r.Sign() <= 0 || r.Cmp(btcec.S256().N) >= 0 {
		return false
	}
	if s.Sign() <= 0 || s.Cmp(btcec.S256().N) >= 0 {
		return false
	}

	hash := sha256.Sum256([]byte(message))
	z := new(big.Int).SetBytes(hash[:])

	w := new(big.Int).ModInverse(s, btcec.S256().N)
	if w == nil {
		return false
	}

	u1 := new(big.Int).Mul(z, w)
	u1.Mod(u1, btcec.S256().N)
	u2 := new(big.Int).Mul(r, w)
	u2.Mod(u2, btcec.S256().N)

	x1, y1 := c.ScalarBaseMult(u1.Bytes())
	x2, y2 := c.ScalarMult(pubkey.X, pubkey.Y, u2.Bytes())
	x, y := c.Add(x1, y1, x2, y2)

	if x.Sign() == 0 && y.Sign() == 0 {
		return false
	}

	v := new(big.Int).Mod(x, btcec.S256().N)
	return v.Cmp(r) == 0
}

func main() {

	var privKeys []*btcec.PrivateKey
	var pubKeys []*btcec.PublicKey
	var Bx, By []*big.Int

	for i := 0; i < btcec.Imax; i++ {
		privKey, err := btcec.NewPrivateKey(btcec.S256())
		if err != nil {
			log.Fatalf("Failed to generate private key: %v", err)
		}
		pubKey := privKey.PubKey()
		privKeys = append(privKeys, privKey)
		pubKeys = append(pubKeys, pubKey)
		Bx = append(Bx, pubKey.X)
		By = append(By, pubKey.Y)
	}

	message := "Hello, ECDSA!"

	start1 := time.Now()
	rs, ss, err := Signbatch(privKeys, message)
	cost1 := time.Since(start1)
	if err != nil {
		log.Fatalf("Failed to sign message: %v", err)
	}
	fmt.Printf("BatchSignature: (rs[0]: %s, ss:[0] %s)\n", rs[0].String(), ss[0].String())

	isValid := VerifySignature(pubKeys[0], message, rs[0], ss[0])
	fmt.Printf("BatchSignature[0] valid: %v\n", isValid)

	start2 := time.Now()
	isValid = VerifySignatureBatch(Bx, By, message, rs, ss)
	cost2 := time.Since(start2)
	fmt.Printf("BatchSignature valid: %v\n", isValid)

	for i := 0; i < btcec.Imax; i++ {

		start3 := time.Now()
		r, s, err := Sign(privKeys[i], message)
		cost3 := time.Since(start3)
		totalcost3 += cost3
		if err != nil {
			log.Fatalf("Failed to sign message: %v", err)
		}
		if i == 0 {
			fmt.Printf("Signature: (r: %s, s: %s)\n", r.String(), s.String())
		}

		start4 := time.Now()
		isValid = VerifySignature(pubKeys[i], message, r, s)
		cost4 := time.Since(start4)
		totalcost4 += cost4
		if i == 0 {
			fmt.Printf("Signature valid: %v\n", isValid)
		}

	}
	fmt.Printf("Sign %d times  cost=[%s]\n", btcec.Imax, totalcost3)
	fmt.Printf("Signbatch %d times  cost=[%s]\n", btcec.Imax, cost1)
	fmt.Printf("%d times Signbatch faster than %f time Sign\n", btcec.Imax, (float64(totalcost3)-float64(cost1))/float64(totalcost3))
	fmt.Printf("VerifySignature %d times  cost=[%s]\n", btcec.Imax, totalcost4)
	fmt.Printf("VerifySignatureBatch %d times  cost=[%s]\n", btcec.Imax, cost2)
	fmt.Printf("%d times VerifySignatureBatch faster than %f time VerifySignature\n", btcec.Imax, (float64(totalcost4)-float64(cost2))/float64(totalcost4))

	fmt.Printf("%d times Signbatch is %f time Sign\n", btcec.Imax, float64(cost1)/float64(totalcost3))
	fmt.Printf("%d times VerifySignatureBatch is %f time VerifySignature\n", btcec.Imax, float64(cost2)/float64(totalcost4))
	fmt.Printf("%d times Signbatch and VerifySignatureBatch is Sign and VerifySignature %f time\n", btcec.Imax, (float64(cost1)+float64(cost2))/(float64(totalcost3)+float64(totalcost4)))
	fmt.Printf("%d times Signbatch and VerifySignatureBatch faster than Sign and VerifySignature %f time\n", btcec.Imax, (float64(totalcost3)+float64(totalcost4)-float64(cost1)-float64(cost2))/(float64(totalcost3)+float64(totalcost4)))

}
