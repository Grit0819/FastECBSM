package main

import (
	"crypto/rand"
	"crypto/sha256"
	"fmt"
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

func Signbatch(priv *btcec.PrivateKey, message string) (rs, ss []*big.Int, err error) {
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
		// Generate a random k
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
			s := new(big.Int).Mul(r, priv.D)
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

func VerifySignatureBatch(pubkey *btcec.PublicKey, message string, rs, ss []*big.Int) bool {

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
	x2s, y2s := c.ScalarMultBatchFieldValSingleQ(pubkey.X, pubkey.Y, u2sBytes)

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

func hashToCurvePoint(input string) (*big.Int, *big.Int, error) {
	var x, y *big.Int

	for {
		h := sha256.New()
		h.Write([]byte(input))
		hash := h.Sum(nil)
		x = new(big.Int).SetBytes(hash)
		x.Mod(x, c.Params().P)

		ySquared := new(big.Int).Exp(x, big.NewInt(3), c.Params().P)
		ySquared.Add(ySquared, c.Params().B)
		ySquared.Mod(ySquared, c.Params().P)

		y = new(big.Int).ModSqrt(ySquared, c.Params().P)
		if y != nil {
			return x, y, nil
		}

		input = fmt.Sprintf("%x", hash)
	}
}

func generateRandomNumber(length int) string {
	const digits = "0123456789"
	b := make([]byte, length)
	for i := range b {
		num, err := rand.Int(rand.Reader, big.NewInt(int64(len(digits))))
		if err != nil {
			panic(err)
		}
		b[i] = digits[num.Int64()]
	}
	return string(b)
}

func generateIDNumber() string {

	addressCode := generateRandomNumber(6)

	now := time.Now()
	startYear := now.Year() - 100
	endYear := now.Year()
	year := startYear + int(generateRandomInt(int64(endYear-startYear+1)))
	month := 1 + int(generateRandomInt(12))
	day := 1 + int(generateRandomInt(28))

	dateCode := fmt.Sprintf("%04d%02d%02d", year, month, day)

	seqCode := generateRandomNumber(3)

	checkCode := generateRandomNumber(1)
	return addressCode + dateCode + seqCode + checkCode
}

func generateRandomInt(max int64) int64 {
	n, err := rand.Int(rand.Reader, big.NewInt(max))
	if err != nil {
		panic(err)
	}
	return n.Int64()
}

func generateTestData(count int, overlapRatio float64) ([]string, []string) {
	data1 := make([]string, count)
	data2 := make([]string, count)

	overlapCount := int(float64(count) * overlapRatio)
	if overlapCount == 0 {
		overlapCount = 1
	}

	usedStrings := make(map[string]bool)

	for i := 0; i < overlapCount; i++ {
		str := generateIDNumber()
		for usedStrings[str] {
			str = generateIDNumber()
		}
		usedStrings[str] = true
		data1[i] = str
		data2[i] = str
	}

	for i := overlapCount; i < count; i++ {
		str := generateIDNumber()
		for usedStrings[str] {
			str = generateIDNumber()
		}
		usedStrings[str] = true
		data1[i] = str
	}

	for i := overlapCount; i < count; i++ {
		str := generateIDNumber()
		for usedStrings[str] {
			str = generateIDNumber()
		}
		usedStrings[str] = true
		data2[i] = str
	}

	for i := 0; i < count; i++ {
		fmt.Printf("data1[%d]: %s, data2[%d]: %s\n", i, data1[i], i, data2[i])
	}

	return data1, data2
}

func countIntersection(set1, set2 []string) int {
	set1Map := make(map[string]struct{})
	for _, item := range set1 {
		set1Map[item] = struct{}{}
	}

	intersectionCount := 0
	for _, item := range set2 {
		if _, found := set1Map[item]; found {
			intersectionCount++
		}
	}

	return intersectionCount
}

func hashPoint(x, y *big.Int) string {
	h := sha256.New()
	h.Write(x.Bytes())
	h.Write(y.Bytes())
	return fmt.Sprintf("%x", h.Sum(nil))
}
func main() {

	var BatchHashToCurvePointTime time.Duration = 0

	const messageCount = 1048576
	const overlapRatio = 0.3

	aliceSet, bobSet := generateTestData(messageCount, overlapRatio)

	intersectionCount := countIntersection(aliceSet, bobSet)

	fmt.Printf("Intersection count between Alice and Bob's messages: %d\n", intersectionCount)

	alicePrivateKey, _ := btcec.NewPrivateKey(c)

	var hashAliceXs []*big.Int
	var hashAliceYs []*big.Int
	for _, x := range aliceSet {
		starthashToCurvePointTime1 := time.Now()
		hashX, hashY, err := hashToCurvePoint(x)
		if err != nil {
			fmt.Println("散列到曲线点时出错:", err)
			return
		}
		costhashToCurvePointTime1 := time.Since(starthashToCurvePointTime1)
		BatchHashToCurvePointTime += costhashToCurvePointTime1
		hashAliceXs = append(hashAliceXs, hashX)
		hashAliceYs = append(hashAliceYs, hashY)

	}

	for i := 0; i < len(hashAliceXs); i++ {
		c.ScalarMult(hashAliceXs[i], hashAliceYs[i], alicePrivateKey.D.Bytes())

	}

	c.ScalarMultBatchFieldVal(hashAliceXs, hashAliceYs, alicePrivateKey.D.Bytes())

}
