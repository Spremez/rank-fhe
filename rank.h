#pragma once
#include "openfhe.h"
#include "encryption.h"
#include <vector>
#include "rotation.h"

// 删除此行LOG2宏定义
// #define LOG2(x) static_cast<int>(std::log2(x))

// 实现矩阵行掩码函数
lbcrypto::Ciphertext<lbcrypto::DCRTPoly> maskRow(
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> c,
    const size_t matrixSize,
    const size_t rowIndex,
    lbcrypto::CryptoContext<lbcrypto::DCRTPoly> cc
);

// 矩阵行求和函数（从sorting-fhe提取）
lbcrypto::Ciphertext<lbcrypto::DCRTPoly> sumRows(
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> c,
    const size_t matrixSize,
    bool maskOutput = false,
    const size_t outputRow = 0,
    // 在第22行添加参数默认值
    lbcrypto::CryptoContext<lbcrypto::DCRTPoly> cc = nullptr
);
Ciphertext<DCRTPoly> innerProduct(
    CryptoContext<DCRTPoly>& cc,
    RotationComposer<10>& rot,
    const Ciphertext<DCRTPoly>& a,
    const Ciphertext<DCRTPoly>& b,
    int vectorSize
);