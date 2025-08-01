#include "rank.h"
#include <vector>
#include "rotation.h"

// 实现矩阵行掩码函数
lbcrypto::Ciphertext<lbcrypto::DCRTPoly> maskRow(
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> c,
    const size_t matrixSize,
    const size_t rowIndex,
    lbcrypto::CryptoContext<lbcrypto::DCRTPoly> cc = nullptr
) {
    std::vector<double> mask(matrixSize * matrixSize, 0.0);
    for (size_t i = 0; i < matrixSize; ++i) {
        mask[matrixSize * rowIndex + i] = 1.0;
    }
    // 在第15行将double向量转换为int64_t向量
    std::vector<int64_t> mask_int(mask.begin(), mask.end());
    auto mask_ct = cc->MakePackedPlaintext(mask_int);
    return cc->EvalMult(c, mask_ct);
}

// 实现矩阵行求和函数
lbcrypto::Ciphertext<lbcrypto::DCRTPoly> sumRows(
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> c,
    const size_t matrixSize,
    bool maskOutput,
    const size_t outputRow,
    lbcrypto::CryptoContext<lbcrypto::DCRTPoly> cc
) {
    for (size_t i = 0; i < LOG2(matrixSize); ++i) {
        c = cc->EvalAdd(c, cc->EvalRotate(c, 1 << (LOG2(matrixSize) + i)));
    }
    if (maskOutput) {
        c = maskRow(c, matrixSize, outputRow, cc);
    }
    return c;
}

Ciphertext<DCRTPoly> innerProduct(CryptoContext<DCRTPoly> &cc, RotationComposer<10> &rot, const Ciphertext<DCRTPoly> &a, const Ciphertext<DCRTPoly> &b, int vectorSize) {
    // 元素-wise 乘法
    auto product = cc->EvalMultAndRelinearize(a, b);
    
    // 旋转累加求和（参考sumColumnsHybrid实现）
    for (int i = 0; i < LOG2(vectorSize); ++i) {
        int rotationStep = 1 << i;
        product = cc->EvalAdd(product, rot.rotate(product, rotationStep));
    }
    
    // 提取结果（第一个槽位）
    std::vector<double> mask(vectorSize, 0.0);
    mask[0] = 1.0;
    Plaintext pmsk = cc->MakeCKKSPackedPlaintext(mask, 1, product->GetLevel());
    return cc->EvalMult(product, pmsk);
}