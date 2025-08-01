#include "encryption.h"
#include "rotation.h"
#include "rank.h"
#include "sign.h"
#include <vector>
#include <random>
#include <chrono>  // 用于时间测量
#include <iomanip>  // 添加此行用于std::setprecision


int main() {
    // 1. 初始化加密上下文（支持1000个槽位）
    lbcrypto::CCParams<lbcrypto::CryptoContextCKKSRNS> parameters;
    parameters.SetSecurityLevel(lbcrypto::HEStd_128_classic); // 添加128位安全级别
    parameters.SetBatchSize(32768);  // SIMD批处理大小 required >100*100
    parameters.SetScalingModSize(40);
    parameters.SetMultiplicativeDepth(50); // 乘法深度
    parameters.SetScalingTechnique(FLEXIBLEAUTO);
    
    auto cc = GenCryptoContext(parameters);
    cc->Enable(lbcrypto::PKE);
    cc->Enable(lbcrypto::KEYSWITCH);
    cc->Enable(lbcrypto::LEVELEDSHE);
    cc->Enable(lbcrypto::ADVANCEDSHE);

    // 2. 生成密钥（添加必要旋转步长）
    auto keyPair = cc->KeyGen();
    uint32_t M = 2 * cc->GetRingDimension();
    std::cout << "M: " << M << std::endl;
    // 添加乘法评估密钥生成
    cc->EvalMultKeyGen(keyPair.secretKey);
    // 修复：添加负方向旋转密钥
    std::vector<int> rotIndices = {-1, -2, -4, -8, -16, -32, -64, 127, 127*2, 127*4,127*8, 127*16, 127*32, 127*64, -128, -(128*2), -(128*4), -(128*8), -(128*16), -(128*32), -(128*64)};    
    cc->EvalRotateKeyGen(keyPair.secretKey, rotIndices);
    std::cout << "旋转密钥生成完成" << std::endl;
    // 3. 准备数据
    std::vector<double> queryVector(10);  // 单个10维查询向量
    std::vector<std::vector<double>> database(100, std::vector<double>(10));  // 100个10维向量

    // 添加随机数生成器，确保分量值在(-10, 10)之间
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 20.0);

    // 生成查询向量
for (auto& val : queryVector) val = dist(gen);

// 对查询向量进行正则化
double queryNorm = 0.0;
for (const auto& val : queryVector) queryNorm += val * val;
queryNorm = std::sqrt(2*queryNorm);
for (auto& val : queryVector) val /= queryNorm;

// 填充数据库向量并正则化，为了接下来compare不同向量的内积结果，需要使得内积结果落在[0,1]
for (auto& vec : database) {
    for (auto& val : vec) val = dist(gen);
    
    // 对当前数据库向量进行正则化
    double norm = 0.0;
    for (const auto& val : vec) norm += val * val;
    norm = std::sqrt(2*norm);
    for (auto& val : vec) val /= norm;
}

    // 4. 数据矩阵化 - 每个10维向量占据128个slot（为了容纳100个向量内积），前10个为值，后118个为0
    std::vector<double> flatDatabase;
    for (const auto& vec : database) {
        flatDatabase.insert(flatDatabase.end(), vec.begin(), vec.end());  // 添加前10个slot的值
        flatDatabase.insert(flatDatabase.end(), 118, 0.0);  // 填充剩余118个slot为0
    }
    
    // 5. 查询向量处理 - 每个10维向量占据128个slot，前10个为值，后118个为0
    std::vector<double> repeatedQuery;
    for (int i = 0; i < 100; ++i) {
        repeatedQuery.insert(repeatedQuery.end(), queryVector.begin(), queryVector.end());  // 添加前10个slot的值
        repeatedQuery.insert(repeatedQuery.end(), 118, 0.0);  // 填充剩余118个slot为0
    }
    
    // 添加：计算明文内积作为验证基准
   /*** std::vector<double> plaintextInnerProducts;
    for (const auto& dbVec : database) {
        double dotProduct = 0.0;
        for (int j = 0; j < 10; ++j) {
            dotProduct += queryVector[j] * dbVec[j];
        }
        plaintextInnerProducts.push_back(dotProduct);
    }*/
    

    // 6. 加密数据
    auto enc = std::make_shared<Encryption>(cc, keyPair.publicKey);
    auto dbCiphertext = enc->encryptInput(flatDatabase);      // 100个向量密文
    auto queryCiphertext = enc->encryptInput(repeatedQuery);  // 100个重复向量密文
    std::cout << "查询向量及数据库向量编码完成" << std::endl;
    // 7. 单次密文乘法实现欧氏距离的分量计算
    auto subCiphertext = cc->EvalSub(dbCiphertext, queryCiphertext);
    auto product = cc->EvalMultAndRelinearize(subCiphertext, subCiphertext);

    // 8. 批量旋转累加（处理100个向量组）
    auto start = high_resolution_clock::now();
    auto summed = product;
    auto cprecomp = cc->EvalFastRotationPrecompute(summed);
    for (int i = 0; i < 4; ++i) {  
        int rotationStep = 1 << i;  // 步长: 1,2,4,8
        // 使用负步长实现向右旋转
        summed = cc->EvalAdd(summed, cc->EvalFastRotation(summed, -rotationStep, M, cprecomp));  
        cprecomp = cc->EvalFastRotationPrecompute(summed);
    }
    

    // 9. 提取100个内积结果并执行rank的过程
    std::vector<double> resultMask(101 * 128, 0.0);
    // 内积结果实际存储在每组第10个槽位（索引9）
    for (int i = 0; i < 100; ++i) resultMask[i * 128 + 9] = 1.0;
    // 显式使用CKKS打包明文编码
    Plaintext pmsk = cc->MakeCKKSPackedPlaintext(resultMask);
    auto IPResult = cc->EvalMult(summed, pmsk);  // 明文乘法
    auto end = high_resolution_clock::now();
    duration<double, std::milli> duration = end - start;
    std::cout << "内积计算完成" <<","<<"ciphertext level:"<<IPResult->GetLevel() << std::endl;
    std::cout << "SIMD内积计算操作耗时: " << duration.count() << " 毫秒" << std::endl;

    // 利用rotation操作将每个密文的内积结果复制到相邻的128个slot内: 第10个到第137个slot, 138-265slot,....
    start = high_resolution_clock::now();
    summed = IPResult;
    for (int i = 0; i < 7; ++i) {  
        int rotationStep = 1 << i;  // 步长: 1,2,4,8,16,32,64
        // 使用负步长实现向右旋转
        cprecomp = cc->EvalFastRotationPrecompute(summed);
        summed = cc->EvalAdd(summed, cc->EvalFastRotation(summed, -rotationStep, M, cprecomp));  
    }
    // 提取每组128个slot的100个slot 结果形式为sum_0 .....sum_0, sum_1...sum_1,.....sum_99...sum_99(每个内积有100个重复)
    for (int i = 0; i < 100; ++i) {
        for (int j = 0; j < 100; ++j) {
            resultMask[i * 128 + 9 + j] = 1.0;
        }
    }
    pmsk = cc->MakeCKKSPackedPlaintext(resultMask);
    auto repeatResult = cc->EvalMult(summed, pmsk);
    end = high_resolution_clock::now();
    duration = end - start;
    std::cout << "repeat内积值完成" <<","<<"ciphertext level:"<<repeatResult->GetLevel() << std::endl;
    std::cout << "SIMD重复内积计算操作耗时: " << duration.count() << " 毫秒" << std::endl;
    //利用rotation操作得到连续的100个内积结果放在相邻的100个slot内，即sum_0,sum_1,...,sum_99
    start = high_resolution_clock::now();
    summed = IPResult;
    for (int i = 0; i < 7; ++i) {
        int rotationStep = (1 << i)*127;  // 步长:127,127*2,127*4,127*8,127*16,127*32,127*64
        // 使用正步长实现向左旋转
        cprecomp = cc->EvalFastRotationPrecompute(summed);
        summed = cc->EvalAdd(summed, cc->EvalFastRotation(summed, rotationStep, M, cprecomp));  
    }
    //提取从第10个slot到第109个slot的值 其形式为sum_0,sum_1,...,sum_99
    for (int i = 1; i < 100; ++i) {
        for (int j = 0; j < 100; ++j) {
            resultMask[i * 128 + 9 + j] = 0.0;
        }
    }
    pmsk = cc->MakeCKKSPackedPlaintext(resultMask);
    auto orderResult = cc->EvalMult(summed, pmsk);
    end = high_resolution_clock::now();
    duration = end - start;
    std::cout << "order内积值完成" <<","<<"ciphertext level:"<<orderResult->GetLevel() << std::endl;
    std::cout << "order内积计算操作耗时: " << duration.count() << " 毫秒" << std::endl;
    //将sum_0,sum_1,...,sum_99复制100次到剩下的100个slot内
    start = high_resolution_clock::now();
    summed = orderResult;
    for (int i = 0; i < 7; ++i) {
        int rotationStep = (1 << i)*128;  // 步长:128,128*2,128*4,128*8,128*16,128*32,128*64
        // 使用正步长实现向左旋转
        cprecomp = cc->EvalFastRotationPrecompute(summed);
        summed = cc->EvalAdd(summed, cc->EvalFastRotation(summed, -rotationStep, M, cprecomp));  
    }
    //提取从每组向量的100个内积 其形式为sum_0,sum_1,...,sum_99
    for (int i = 0; i < 100; ++i) {
        for (int j = 0; j < 100; ++j) {
            resultMask[i * 128 + 9 + j] = 1.0;
        }
    }
    pmsk = cc->MakeCKKSPackedPlaintext(resultMask);
    auto repeatorderResult = cc->EvalMult(summed, pmsk);
    end = high_resolution_clock::now();
    duration = end - start;
    std::cout << "repeat-order内积值完成" <<","<<"ciphertext level:"<<repeatorderResult->GetLevel() << std::endl;
    std::cout << "repeat-order内积计算操作耗时: " << duration.count() << " 毫秒" << std::endl;
    //将repeatResult：sum_0 .....sum_0, sum_1...sum_1,.....sum_99...sum_99(每个内积有100个重复)同repeatorderResult：重复的sum_0,sum_1,...,sum_99相减
    auto subResult = cc->EvalSub(repeatResult, repeatorderResult);
    // Parameters for compositeSign
    int dg = 1;
    int df = 2;

    // Apply compositeSign to subResult 
    start = high_resolution_clock::now();
    auto result = compositeSign<4>(subResult, cc,
                                   SignConfig(CompositeSignConfig(4, dg, df)));
    end = high_resolution_clock::now();
    duration = end - start;
    std::cout << "compare(compositeSign)完成" <<","<<"ciphertext level:"<<result->GetLevel() << std::endl;
    std::cout << "compare(compositeSign)计算操作耗时: " << duration.count() << " 毫秒" << std::endl;
    //将每组100个内积的比较值相加得到每个元素的rank值
    start = high_resolution_clock::now();
    summed =result;
    for (int i = 0; i < 7; ++i) {  
        int rotationStep = 1 << i;  // 步长: 1,2,4,8,16,32,64
        // 使用负步长实现向右旋转
        cprecomp = cc->EvalFastRotationPrecompute(summed);
        summed = cc->EvalAdd(summed, cc->EvalFastRotation(summed, -rotationStep, M, cprecomp));  
    }
    //每个元素的rank值存储在每组的第109个slot内
    //rank值的范围为[-99,99]
    //为了得到rank最小的10个元素，需要将rank值取反，并将所有rank值-90即可得到rank最大的10个元素为非负值其余为负值
    for (int i = 0; i < 100; ++i) {
        resultMask[i * 128 + 108] = -1.0;
        for (int j = 0; j < 99; ++j) {
            resultMask[i * 128 + 9 + j] = 0.0;
        }
    }
    pmsk = cc->MakeCKKSPackedPlaintext(resultMask);
    auto rankResult = cc->EvalMult(summed, pmsk);
    end = high_resolution_clock::now();
    duration = end - start;
    std::cout << "rank计算完成" <<","<<"ciphertext level:"<<summed->GetLevel() << std::endl;
    std::cout << "rank计算操作耗时: " << duration.count() << " 毫秒" << std::endl;
    //将rankresult对应slot-90
    for (int i = 0; i < 100; ++i) {
        resultMask[i * 128 + 108] = -89.0;
    }
    pmsk = cc->MakeCKKSPackedPlaintext(resultMask);
    rankResult = cc->EvalAdd(rankResult, pmsk);

    // Apply compositeSign to rankResult 
    result = compositeSign<4>(rankResult, cc,
                                   SignConfig(CompositeSignConfig(4, dg, df)));
    std::cout << "TOP10计算完成" <<","<<"ciphertext level:"<<result->GetLevel() << std::endl;

    // 10. 解密结果
    Plaintext plainResult;
    //cc->Decrypt(keyPair.secretKey, IPResult, &plainResult);
    //std::vector<double> decrypted = plainResult->GetRealPackedValue();
    auto debugEnc = std::make_shared<DebugEncryption>(cc, keyPair);
    
    /*** std::cout << "内积计算结果：\n";
    // 添加：设置输出精度
    std::cout << std::fixed << std::setprecision(4);
   
    for (int i = 0; i < 100; ++i) {
        // 添加：计算绝对误差
        double error = std::abs(decrypted[i*128+9] - plaintextInnerProducts[i]);
        std::cout << "第" << i << "个内积: "
                  << "同态结果=" << decrypted[i*128+9] << ", "
                  << "明文结果=" << plaintextInnerProducts[i] << ", "
                  << "误差=" << error << ", "
                  << (error < 1e-3 ? "验证通过" : "验证失败") << "\n";
    }
    ***/

    cc->Decrypt(keyPair.secretKey, result, &plainResult);
    std::vector<double> decrypted = plainResult->GetRealPackedValue();
    std::cout << "内积计算结果：\n";
    // 添加：设置输出精度
    std::cout << std::fixed << std::setprecision(2);
    // 统计非负值个数并记录对应的i值
    int nonNegativeCount = 0;
    std::vector<int> nonNegativeIndices;
    
    for (int i = 0; i < 100; ++i) {
        double value = decrypted[i*128+108];
        if (value >= 0) {
            nonNegativeCount++;
            nonNegativeIndices.push_back(i);
        }
    }
    
    // 输出结果
    std::cout << "非负值个数: " << nonNegativeCount << std::endl;
    std::cout << "非负值对应的i值: ";
    for (size_t j = 0; j < nonNegativeIndices.size(); ++j) {
        std::cout << nonNegativeIndices[j];
        if (j < nonNegativeIndices.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << std::endl;

    return 0;
}