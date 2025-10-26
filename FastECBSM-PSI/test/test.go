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
		go c.ScalarMult(c1xs[i], c1ys[i], privKey.D.Bytes()) // go并行 解密
	}

	btcec.WG.Wait()
	cost10 := time.Since(startW)
	costW = cost10

	for i := 0; i < btcec.TestNum; i++ {
		cipher1 := btcec.Encrypt(pubKey, messages[btcec.TestNum-1-i]) //使用公钥对消息进行同态加密

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
	var scalarTime time.Duration = 0
	var BatchScalarTime time.Duration = 0
	var hashToCurvePointTime time.Duration = 0
	var BatchHashToCurvePointTime time.Duration = 0

	const messageCount = 1048576
	const overlapRatio = 0.3

	aliceSet, bobSet := generateTestData(messageCount, overlapRatio)

	intersectionCount := countIntersection(aliceSet, bobSet)

	fmt.Printf("Intersection count between Alice and Bob's messages: %d\n", intersectionCount)

	alicePrivateKey, _ := btcec.NewPrivateKey(c)
	bobPrivateKey, _ := btcec.NewPrivateKey(c)

	start1 := time.Now()

	var AliceXs []*big.Int
	var AliceYs []*big.Int
	for _, x := range aliceSet {
		starthashToCurvePointTime1 := time.Now()
		hashX, hashY, err := hashToCurvePoint(x)
		if err != nil {
			fmt.Println("散列到曲线点时出错:", err)
			return
		}
		costhashToCurvePointTime1 := time.Since(starthashToCurvePointTime1)
		hashToCurvePointTime += costhashToCurvePointTime1
		startScalarTime1 := time.Now()
		ax, ay := c.ScalarMult(hashX, hashY, alicePrivateKey.D.Bytes())
		costScalarTime1 := time.Since(startScalarTime1)
		scalarTime += costScalarTime1
		AliceXs = append(AliceXs, ax)
		AliceYs = append(AliceYs, ay)

	}

	var BobXs []*big.Int
	var BobYs []*big.Int
	for _, y := range bobSet {
		starthashToCurvePointTime2 := time.Now()
		hashX, hashY, err := hashToCurvePoint(y)
		if err != nil {
			fmt.Println("散列到曲线点时出错:", err)
			return
		}
		costhashToCurvePointTime2 := time.Since(starthashToCurvePointTime2)
		hashToCurvePointTime += costhashToCurvePointTime2
		startScalarTime2 := time.Now()
		bx, by := c.ScalarMult(hashX, hashY, bobPrivateKey.D.Bytes())
		costScalarTime2 := time.Since(startScalarTime2)
		scalarTime += costScalarTime2
		BobXs = append(BobXs, bx)
		BobYs = append(BobYs, by)

	}

	var bobSharedSecretXs []*big.Int
	var bobSharedSecretYs []*big.Int
	for i := 0; i < len(AliceXs); i++ {
		startScalarTime3 := time.Now()
		sharedX, sharedY := c.ScalarMult(AliceXs[i], AliceYs[i], bobPrivateKey.D.Bytes())
		costScalarTime3 := time.Since(startScalarTime3)
		scalarTime += costScalarTime3
		bobSharedSecretXs = append(bobSharedSecretXs, sharedX)
		bobSharedSecretYs = append(bobSharedSecretYs, sharedY)
	}

	inverseSa := new(big.Int).ModInverse(alicePrivateKey.D, c.Params().N)

	var aliceSharedSecretXs []*big.Int
	var aliceSharedSecretYs []*big.Int
	for i := 0; i < len(bobSharedSecretXs); i++ {
		startScalarTime4 := time.Now()
		sharedX, sharedY := c.ScalarMult(bobSharedSecretXs[i], bobSharedSecretYs[i], inverseSa.Bytes())
		costScalarTime4 := time.Since(startScalarTime4)
		scalarTime += costScalarTime4
		aliceSharedSecretXs = append(aliceSharedSecretXs, sharedX)
		aliceSharedSecretYs = append(aliceSharedSecretYs, sharedY)
	}

	cost3 := time.Since(start1)

	res1 := 0

	startcmp1 := time.Now()
	sharedSet := make(map[string]struct{})

	for i := 0; i < len(BobXs); i++ {
		sharedSet[BobXs[i].String()] = struct{}{}
	}

	for i := 0; i < len(aliceSharedSecretXs); i++ {
		if _, found := sharedSet[aliceSharedSecretXs[i].String()]; found {
			res1++
		}
	}
	costcmp1 := time.Since(startcmp1)
	cost1 := time.Since(start1)

	start2 := time.Now()

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
	BatchStartScalarTime1 := time.Now()
	alicePublicPointsBatchX, alicePublicPointsBatchY := c.ScalarMultBatchFieldVal(hashAliceXs, hashAliceYs, alicePrivateKey.D.Bytes())
	BatchCostScalarTime1 := time.Since(BatchStartScalarTime1)
	BatchScalarTime += BatchCostScalarTime1

	var hashBobXs []*big.Int
	var hashBobYs []*big.Int
	for _, y := range bobSet {
		starthashToCurvePointTime2 := time.Now()
		hashBobX, hashBobY, err := hashToCurvePoint(y)
		if err != nil {
			fmt.Println("散列到曲线点时出错:", err)
			return
		}
		costhashToCurvePointTime2 := time.Since(starthashToCurvePointTime2)
		BatchHashToCurvePointTime += costhashToCurvePointTime2
		hashBobXs = append(hashBobXs, hashBobX)
		hashBobYs = append(hashBobYs, hashBobY)

	}
	BatchStartScalarTime2 := time.Now()
	bobPublicPointsBatchX, _ := c.ScalarMultBatchFieldVal(hashBobXs, hashBobYs, bobPrivateKey.D.Bytes())
	BatchCostScalarTime2 := time.Since(BatchStartScalarTime2)
	BatchScalarTime += BatchCostScalarTime2

	BatchStartScalarTime3 := time.Now()
	BatchbobSharedSecretsX, BatchbobSharedSecretsY := c.ScalarMultBatchFieldVal(alicePublicPointsBatchX, alicePublicPointsBatchY, bobPrivateKey.D.Bytes())
	BatchCostScalarTime3 := time.Since(BatchStartScalarTime3)
	BatchScalarTime += BatchCostScalarTime3

	inverseSaBatch := new(big.Int).ModInverse(alicePrivateKey.D, c.Params().N)

	BatchStartScalarTime4 := time.Now()
	BatchaliceSharedSecretsX, _ := c.ScalarMultBatchFieldVal(BatchbobSharedSecretsX, BatchbobSharedSecretsY, inverseSaBatch.Bytes())
	BatchCostScalarTime4 := time.Since(BatchStartScalarTime4)
	BatchScalarTime += BatchCostScalarTime4

	cost4 := time.Since(start2)

	res := 0
	startcmp2 := time.Now()
	sharedSetBatch := make(map[string]struct{})

	for i := 0; i < len(bobPublicPointsBatchX); i++ {
		sharedSetBatch[bobPublicPointsBatchX[i].String()] = struct{}{}
	}

	for i := 0; i < len(BatchaliceSharedSecretsX); i++ {
		if _, found := sharedSetBatch[BatchaliceSharedSecretsX[i].String()]; found {
			res++
		}
	}
	costcmp2 := time.Since(startcmp2)

	cost2 := time.Since(start2)

	fmt.Printf("singal Intersection count between Ga x Sb and Gb × Sb : %d\n", res1)
	fmt.Printf("Batch Intersection count between Ga x Sb and Gb × Sb : %d\n", res)

	//时间消耗
	fmt.Printf("PSIsingal's HashToPoints cost=[%s]\n", hashToCurvePointTime)
	fmt.Printf("PSIsingal's cmp cost=[%s]\n", costcmp1)
	fmt.Printf("PSIsingal's TotalScalarTime cost=[%s]\n", scalarTime)
	fmt.Printf("PSIsingal's HashToPoints and TotalScalarTime cost=[%s]\n", cost3)
	fmt.Printf("PSIsingal %d times TotalScalarTime and cmp cost=[%s]\n", messageCount, scalarTime+costcmp1)
	fmt.Printf("PSIsingal %d times  cost=[%s]\n", messageCount, cost1)
	fmt.Printf("\n")
	fmt.Printf("PSIbatch's HashToPoints cost=[%s]\n", BatchHashToCurvePointTime)
	fmt.Printf("PSIbatch's cmp cost=[%s]\n", costcmp2)
	fmt.Printf("PSIsingal's TotalScalarTime cost=[%s]\n", BatchScalarTime)
	fmt.Printf("PSIbatch's HashToPoints and TotalScalarTime cost=[%s]\n", cost4)
	fmt.Printf("PSIbatch %d times TotalScalarTime and cmp cost=[%s]\n", messageCount, BatchScalarTime+costcmp2)
	fmt.Printf("PSIbatch %d times  cost=[%s]\n", messageCount, cost2)

	fmt.Printf("\n")

	fmt.Printf("%d times TotalScalarTime  PSIbatch faster %f time  than PSIsingal\n", messageCount, float64(scalarTime-BatchScalarTime)/float64(scalarTime))
	fmt.Printf("%d times TotalScalarTime PSIbatch is %f time PSIsingal\n", messageCount, float64(BatchScalarTime)/float64(scalarTime))

	fmt.Printf("\n")

	fmt.Printf("%d times ExceptHashtoPoints PSIbatch faster %f time  than PSIsingal\n", messageCount, float64(scalarTime+costcmp1-BatchScalarTime-costcmp2)/float64(scalarTime+costcmp1))
	fmt.Printf("%d times ExceptHashtoPoints PSIbatch is %f time PSIsingal\n", messageCount, float64(BatchScalarTime+costcmp2)/float64(scalarTime+costcmp1))

	fmt.Printf("\n")

	fmt.Printf("%d times PSIbatch faster %f time  than PSIsingal\n", messageCount, float64(cost1-cost2)/float64(cost1))
	fmt.Printf("%d times PSIbatch is %f time PSIsingal\n", messageCount, float64(cost2)/float64(cost1))

}
