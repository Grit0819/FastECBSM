package main

import (
	"crypto/rand"
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

func main() {

	testnum := btcec.TestNum
	var msgmax, _ = new(big.Int).SetString("115792089237316195423570985008687907853269984665640564039457584007913129639935", 10)
	var messages [btcec.TestNum]*big.Int
	k := make([][]byte, len(messages))
	for i := 0; i < testnum; i++ {
		msg, _ := rand.Int(rand.Reader, msgmax)
		if i%2 == 0 {
			messages[i] = msg
		} else {
			messages[i] = msg.Neg(msg)
		}
		k[i] = messages[i].Bytes()
	}

	c := btcec.S256()

	rx := make([]*big.Int, btcec.Imax, btcec.Imax)
	ry := make([]*big.Int, btcec.Imax, btcec.Imax)

	for i := 0; i < btcec.Imax; i++ {
		rx[i], ry[i] = c.ScalarBaseMult(messages[i].Bytes())

	}

	var resx1, resy1 *big.Int
	var resx5, resy5 *big.Int

	start1 := time.Now()
	for i := 0; i < btcec.Imax; i++ {
		resx1, resy1 = c.ScalarMult(rx[i], ry[i], k[i])
	}
	cost1 := time.Since(start1)

	start3 := time.Now()
	resx3, resy3 := c.ScalarMultBatchFieldValBatchQ(rx, ry, k)

	cost3 := time.Since(start3)

	start4 := time.Now()
	resx4, resy4 := c.ScalarBaseMultBatch(k)

	cost4 := time.Since(start4)

	start5 := time.Now()
	for i := 0; i < btcec.Imax; i++ {
		resx5, resy5 = c.ScalarBaseMult(k[i])
	}
	cost5 := time.Since(start5)

	fmt.Println(resx3[btcec.Imax-1], resy3[btcec.Imax-1])
	fmt.Println(resx1, resy1)
	fmt.Println(resx4[btcec.Imax-1], resy4[btcec.Imax-1])
	fmt.Println(resx5, resy5)

	fmt.Printf("-----------------------------------\n")
	fmt.Printf("ScalarMult %d times  cost=[%s]\n", btcec.Imax, cost1)
	fmt.Printf("ScalarMultBatchFieldVal %d times  cost=[%s]\n", btcec.Imax, cost3)
	fmt.Printf("ScalarMultBatchFieldVal %d times is ScalarMult %f\n", btcec.Imax, float64(cost3)/float64(cost1))
	fmt.Printf("ScalarMultBatchFieldVal %d times faster than ScalarMult %f\n", btcec.Imax, float64(cost1-cost3)/float64(cost1))
	fmt.Printf("-----------------------------------\n")
	fmt.Printf("ScalarBaseMult %d times  cost=[%s]\n", btcec.Imax, cost5)
	fmt.Printf("ScalarBaseMultBatch %d times  cost=[%s]\n", btcec.Imax, cost4)
	fmt.Printf("ScalarBaseMultBatch %d times is ScalarBaseMult %f\n", btcec.Imax, float64(cost4)/float64(cost5))
	fmt.Printf("ScalarBaseMultBatch %d times faster than ScalarBaseMult %f\n", btcec.Imax, float64(cost5-cost4)/float64(cost5))
	fmt.Printf("-----------------------------------\n")

}
