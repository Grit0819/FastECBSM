package btcec

import (
	"compress/zlib"
	"encoding/base64"
	"encoding/binary"
	"io/ioutil"
	"strings"
)

//go:generate go run -tags gensecp256k1 genprecomps.go

// loadS256BytePoints decompresses and deserializes the pre-computed byte points
// used to accelerate scalar base multiplication for the secp256k1 curve.  This
// approach is used since it allows the compile to use significantly less ram
// and be performed much faster than it is with hard-coding the final in-memory
// data structure.  At the same time, it is quite fast to generate the in-memory
// data structure at init time with this approach versus computing the table.
func loadS256BytePoints() error {
	// There will be no byte points to load when generating them.
	//secp256k1BytePoints 中获取的预先计算的字节点数据
	bp := secp256k1BytePoints
	if len(bp) == 0 {
		return nil
	}
	//解压用于加速标量基点乘法的预先计算表。
	// Decompress the pre-computed table used to accelerate scalar base
	// multiplication.
	//解码器使用了base64.StdEncoding进行解码，读取bp中的点
	decoder := base64.NewDecoder(base64.StdEncoding, strings.NewReader(bp))
	//解压
	r, err := zlib.NewReader(decoder)
	if err != nil {
		return err
	}
	//读取解压的数据
	serialized, err := ioutil.ReadAll(r)
	if err != nil {
		return err
	}

	// Deserialize the precomputed byte points and set the curve to them.
	//反序列化预先计算的字节点，并将曲线设置为它们。
	//设置最终的预计算的点 bytePoints

	offset := 0
	//定义一个包含 32 个字节的三维数组，用于存储字节点的坐标值
	//32个字节来表示一个坐标点，256表示每个窗口都有256个字节点，3表示坐标为x,y,z
	var bytePoints [32][256][3]FieldVal
	//循环处理每一个字节块
	for byteNum := 0; byteNum < 32; byteNum++ {
		// All points in this window.
		for i := 0; i < 256; i++ {
			px := &bytePoints[byteNum][i][0]
			py := &bytePoints[byteNum][i][1]
			pz := &bytePoints[byteNum][i][2]
			//接下来的三个循环，每个循环将从反序列化后的数据中读取 4 个字节，并将其作为无符号整数存储到当前字节点的坐标值中
			for i := 0; i < 10; i++ {
				px.n[i] = binary.LittleEndian.Uint32(serialized[offset:])
				offset += 4
			}
			for i := 0; i < 10; i++ {
				py.n[i] = binary.LittleEndian.Uint32(serialized[offset:])
				offset += 4
			}
			for i := 0; i < 10; i++ {
				pz.n[i] = binary.LittleEndian.Uint32(serialized[offset:])
				offset += 4
			}
		}
	}
	secp256k1.bytePoints = &bytePoints
	return nil
}
