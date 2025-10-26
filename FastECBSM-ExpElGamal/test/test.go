package main

import (
	"crypto/rand"
	"fmt"
	"log"
	"math/big"

	"time"

	btcec "github.com/cuckoobtcec"
)

// 时间变量的申明
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

var curve = btcec.S256()

type Cipher struct {
	C1x, C1y *big.Int
	C2x, C2y *big.Int
}

func ExpElGamalEncrypt(pubKey *btcec.PublicKey, message *big.Int) (*Cipher, error) {
	r, err := rand.Int(rand.Reader, curve.N)
	if err != nil {
		return nil, fmt.Errorf("failed to generate random r: %v", err)
	}

	C1x, C1y := curve.ScalarBaseMult(r.Bytes())
	hx, hy := curve.ScalarMult(pubKey.X, pubKey.Y, r.Bytes())
	gmx, gmy := curve.ScalarBaseMult(message.Bytes())
	C2x, C2y := curve.Add(gmx, gmy, hx, hy)

	return &Cipher{
		C1x: C1x,
		C1y: C1y,
		C2x: C2x,
		C2y: C2y,
	}, nil
}

func ExpElGamalDecrypt(privKey *btcec.PrivateKey, cipher *Cipher) (*big.Int, error) {
	skC1x, skC1y := curve.ScalarMult(cipher.C1x, cipher.C1y, privKey.D.Bytes())
	invC1y := new(big.Int)
	invC1y.Add(curve.P, invC1y)
	invC1y.Sub(invC1y, skC1y)

	gmx, _ := curve.Add(cipher.C2x, cipher.C2y, skC1x, invC1y)

	message := new(big.Int).SetBytes(gmx.Bytes())
	return message, nil
}

func ExpElGamalEncryptBatch(pubKey *btcec.PublicKey, messages []*big.Int) ([]*Cipher, error) {
	n := len(messages)
	ciphers := make([]*Cipher, n)
	rs := make([]*big.Int, n)
	rsBytes := make([][]byte, n)
	messagesBytes := make([][]byte, n)

	for i := 0; i < n; i++ {
		r, err := rand.Int(rand.Reader, curve.N)
		if err != nil {
			return nil, fmt.Errorf("failed to generate random r: %v", err)
		}
		rs[i] = r
		rsBytes[i] = r.Bytes()
	}

	C1xs, C1ys := curve.ScalarBaseMultBatchFieldVal(rsBytes)

	hxs, hys := curve.ScalarMultBatchFieldValBatchK(pubKey.X, pubKey.Y, rsBytes)
	gmxs := make([]*big.Int, n)
	gmys := make([]*big.Int, n)
	for i := 0; i < n; i++ {
		messagesBytes[i] = messages[i].Bytes() // g^m
	}
	gmxs, gmys = curve.ScalarBaseMultBatchFieldVal(messagesBytes)

	for i := 0; i < n; i++ {
		C2x, C2y := curve.Add(gmxs[i], gmys[i], hxs[i], hys[i])
		ciphers[i] = &Cipher{
			C1x: C1xs[i],
			C1y: C1ys[i],
			C2x: C2x,
			C2y: C2y,
		}
	}

	return ciphers, nil
}

func ExpElGamalDecryptBatch(privKey *btcec.PrivateKey, ciphers []*Cipher) ([]*big.Int, error) {
	n := len(ciphers)
	messages := make([]*big.Int, n)
	C1xs := make([]*big.Int, n)
	C1ys := make([]*big.Int, n)

	for i := 0; i < n; i++ {
		C1xs[i] = ciphers[i].C1x
		C1ys[i] = ciphers[i].C1y
	}

	skC1xs, skC1ys := curve.ScalarMultBatchFieldVal(C1xs, C1ys, privKey.D.Bytes())

	for i := 0; i < n; i++ {
		invC1y := new(big.Int)
		invC1y.Add(curve.P, invC1y)
		invC1y.Sub(invC1y, skC1ys[i])

		gmx, _ := curve.Add(ciphers[i].C2x, ciphers[i].C2y, skC1xs[i], invC1y)
		messages[i] = new(big.Int).SetBytes(gmx.Bytes())
	}

	return messages, nil
}

func main() {

	privKey, err := btcec.NewPrivateKey(curve)
	if err != nil {
		log.Fatalf("Failed to generate private key: %v", err)
	}
	pubKey := privKey.PubKey()

	numMessages := btcec.TestNum
	messages := make([]*big.Int, numMessages)
	for i := 0; i < numMessages; i++ {
		messages[i] = big.NewInt(int64(i + 1))
	}
	var decryptedGmx *big.Int

	fmt.Println("Starting single encryption/decryption for %v times...", btcec.TestNum)

	ciphersSingle := make([]*Cipher, numMessages)
	startSingleEncrypt := time.Now()
	for i := 0; i < numMessages; i++ {
		ciphersSingle[i], _ = ExpElGamalEncrypt(pubKey, messages[i])
	}
	singleEncryptDuration := time.Since(startSingleEncrypt)

	startSingleDecrypt := time.Now()
	for i := 0; i < numMessages; i++ {
		decryptedGmx, _ = ExpElGamalDecrypt(privKey, ciphersSingle[i])
	}
	singleDecryptDuration := time.Since(startSingleDecrypt)

	expectedGmx, _ := curve.ScalarBaseMult(messages[btcec.TestNum-1].Bytes())
	if decryptedGmx.Cmp(expectedGmx) != 0 {
		log.Fatalf("Decryption failed for message %d", btcec.TestNum)
	} else {
		fmt.Printf("最后一次解密验证正确\n")
	}

	fmt.Printf("Single encryption time for %d messages: %s\n", btcec.TestNum, singleEncryptDuration)
	fmt.Printf("Single decryption time for %d messages: %s\n", btcec.TestNum, singleDecryptDuration)

	fmt.Println("Starting batch encryption/decryption for times : ", btcec.TestNum)

	startBatchEncrypt := time.Now()
	ciphersBatch, _ := ExpElGamalEncryptBatch(pubKey, messages)
	batchEncryptDuration := time.Since(startBatchEncrypt)

	startBatchDecrypt := time.Now()
	decryptedBatchGmxs, _ := ExpElGamalDecryptBatch(privKey, ciphersBatch)
	batchDecryptDuration := time.Since(startBatchDecrypt)

	if decryptedBatchGmxs[btcec.TestNum-1].Cmp(expectedGmx) != 0 {
		log.Fatalf("Decryption failed for message %d", btcec.TestNum)
	} else {
		fmt.Printf("最后一次解密验证正确\n")
	}
	fmt.Printf("Batch encryption time for %d messages: %s\n", btcec.TestNum, batchEncryptDuration)
	fmt.Printf("Batch decryption time for %d messages: %s\n", btcec.TestNum, batchDecryptDuration)

	fmt.Println("\n------ Time comparison between single and batch operations ------")

	encryptSpeedup := (float64(singleEncryptDuration) - float64(batchEncryptDuration)) / float64(singleEncryptDuration)
	decryptSpeedup := (float64(singleDecryptDuration) - float64(batchDecryptDuration)) / float64(singleDecryptDuration)
	totalSpeedup := (float64(singleDecryptDuration) + float64(singleEncryptDuration) - float64(batchEncryptDuration) - float64(batchDecryptDuration)) / (float64(singleDecryptDuration) + float64(singleEncryptDuration))

	fmt.Printf("Single encryption vs Batch encryption: %s vs %s\n", singleEncryptDuration, batchEncryptDuration)
	fmt.Printf("Single decryption vs Batch decryption: %s vs %s\n", singleDecryptDuration, batchDecryptDuration)
	fmt.Printf("Batch encryption is %f times faster than single encryption\n", encryptSpeedup)
	fmt.Printf("Batch decryption is %f times faster than single decryption\n", decryptSpeedup)

	fmt.Printf("Batch decryption and encryption is %f times faster than single decryption and encryption\n", totalSpeedup)

	fmt.Println("Encryption and decryption validation successful!")

	fmt.Println("\n\nStarting HomoAdd HomoMul BatchHomoAdd BatchHomoMul for times : ", btcec.TestNum)

	ciphers1 := make([]*btcec.Cipher, btcec.TestNum)
	var testHomoAdd *btcec.Cipher
	var testHomoMul *btcec.Cipher
	for i := 0; i < btcec.TestNum; i++ {
		ciphers1[i] = btcec.Encrypt(pubKey, messages[i])
	}

	ciphers2 := make([]*btcec.Cipher, btcec.TestNum)
	for i, j := btcec.TestNum-1, 0; i >= 0; i, j = i-1, j+1 {
		ciphers2[j] = btcec.Encrypt(pubKey, messages[i])
	}

	for i := 0; i < btcec.TestNum; i++ {
		start4 := time.Now()
		testHomoAdd = btcec.HomoAdd(ciphers1[i], ciphers2[i])
		cost4 := time.Since(start4)
		start5 := time.Now()
		btcec.HomoAddField(ciphers1[i], ciphers2[i])
		cost5 := time.Since(start5)
		start6 := time.Now()
		testHomoMul = btcec.HomoMul(ciphers1[i], messages[i])
		cost6 := time.Since(start6)
		ecchaddall = ecchaddall + cost4
		ecchaddfieldall = ecchaddfieldall + cost5
		ecchmulall = ecchmulall + cost6
	}

	start := time.Now()
	testBatchHomoAdd := btcec.BatchHomoAdd(ciphers1, ciphers2)
	batchAddDuration := time.Since(start)

	start = time.Now()
	testBatchHomoMul := btcec.BatchHomoMul(ciphers1, messages)
	batchMulDuration := time.Since(start)

	if testHomoAdd.C1x.Cmp(testBatchHomoAdd[btcec.TestNum-1].C1x) != 0 || testHomoAdd.C1y.Cmp(testBatchHomoAdd[btcec.TestNum-1].C1y) != 0 ||
		testHomoAdd.C2x.Cmp(testBatchHomoAdd[btcec.TestNum-1].C2x) != 0 || testHomoAdd.C2y.Cmp(testBatchHomoAdd[btcec.TestNum-1].C2y) != 0 {
		fmt.Printf("BatchHomoAdd:Mismatch at index %d: test1 and test2 results are not the same\n", btcec.TestNum-1)
	} else {
		fmt.Printf("BatchHomoAdd:Test passed at index %d: test1 and test2 results are the same\n", btcec.TestNum-1)
	}
	if testHomoMul.C1x.Cmp(testBatchHomoMul[btcec.TestNum-1].C1x) != 0 ||
		testHomoMul.C1y.Cmp(testBatchHomoMul[btcec.TestNum-1].C1y) != 0 ||
		testHomoMul.C2x.Cmp(testBatchHomoMul[btcec.TestNum-1].C2x) != 0 ||
		testHomoMul.C2y.Cmp(testBatchHomoMul[btcec.TestNum-1].C2y) != 0 {
		fmt.Printf("BatchHomoMul:Mismatch at index %d: test1 and test2 results are not the same\n", btcec.TestNum-1)
	} else {
		fmt.Printf("BatchHomoMul:Test passed at index %d: test1 and test2 results are the same\n", btcec.TestNum-1)
	}

	addSpeedup := (float64(ecchaddall) - float64(batchAddDuration)) / float64(ecchaddall)
	mulSpeedup := (float64(ecchmulall) - float64(batchMulDuration)) / float64(ecchmulall)

	fmt.Printf("\n\nTest count (btcec.TestNum): %d\n", btcec.TestNum)

	fmt.Printf("Single HomoAdd duration: %s, Single HomoAddField duration: %s, Batch HomoAdd duration: %s\n", ecchaddall, ecchaddfieldall, batchAddDuration)
	fmt.Printf("Single HomoMul duration: %s, Batch HomoMul duration: %s\n", ecchmulall, batchMulDuration)

	fmt.Printf("Batch HomoAdd is %f times faster than single HomoAdd\n", addSpeedup)
	fmt.Printf("Batch HomoMul is %f times faster than single HomoMul\n", mulSpeedup)

}
