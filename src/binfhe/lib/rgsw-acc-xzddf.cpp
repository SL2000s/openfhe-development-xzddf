// TODO: now copy of AP/DM, change to XZDDF

#include "rgsw-acc-xzddf.h"

#include <string>

namespace lbcrypto {
/*
// Key generation as described in Section 4 of https://eprint.iacr.org/2014/816
RingGSWACCKey RingGSWAccumulatorXZDDF::KeyGenAcc(const std::shared_ptr<RingGSWCryptoParams>& params,
                                              const NativePoly& skNTT, ConstLWEPrivateKey& LWEsk) const {
    auto sv{LWEsk->GetElement()};
    auto mod{sv.GetModulus().ConvertToInt<int32_t>()};
    auto modHalf{mod >> 1};
    uint32_t n(sv.GetLength());
    int32_t baseR(params->GetBaseR());
    const auto& digitsR = params->GetDigitsR();
    RingGSWACCKey ek    = std::make_shared<RingGSWACCKeyImpl>(n, baseR, digitsR.size());


#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(n))
    for (uint32_t i = 0; i < n; ++i) {
        for (int32_t j = 1; j < baseR; ++j) {
            for (size_t k = 0; k < digitsR.size(); ++k) {
                auto s{sv[i].ConvertToInt<int32_t>()};
                (*ek)[i][j][k] =
                    KeyGenDM(params, skNTT, (s > modHalf ? s - mod : s) * j * digitsR[k].ConvertToInt<int32_t>());
            }
        }
    }
    return ek;
}


void RingGSWAccumulatorXZDDF::EvalAcc(const std::shared_ptr<RingGSWCryptoParams>& params, ConstRingGSWACCKey& ek,
                                   RLWECiphertext& acc, const NativeVector& a) const {
    NativeInteger baseR{params->GetBaseR()};
    auto q       = params->Getq();
    auto digitsR = params->GetDigitsR().size();
    uint32_t n   = a.GetLength();

    for (uint32_t i = 0; i < n; ++i) {
        auto aI = NativeInteger(0).ModSubFast(a[i], q);
        for (size_t k = 0; k < digitsR; ++k, aI /= baseR) {
            auto a0 = (aI.Mod(baseR)).ConvertToInt<uint32_t>();
            if (a0)
                AddToAccDM(params, (*ek)[i][a0][k], acc);
        }
    }
}

// AP Accumulation as described in https://eprint.iacr.org/2020/086
void RingGSWAccumulatorXZDDF::AddToAccDM(const std::shared_ptr<RingGSWCryptoParams>& params, ConstRingGSWEvalKey& ek,   // TODO: change DM to XZDDF
                                      RLWECiphertext& acc) const {
    std::vector<NativePoly> ct(acc->GetElements());
    ct[0].SetFormat(Format::COEFFICIENT);
    ct[1].SetFormat(Format::COEFFICIENT);

    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG2{(params->GetDigitsG() - 1) << 1};
    std::vector<NativePoly> dct(digitsG2, NativePoly(params->GetPolyParams(), Format::COEFFICIENT, true));

    SignedDigitDecompose(params, ct, dct);

#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(digitsG2))
    for (uint32_t j = 0; j < digitsG2; ++j)
        dct[j].SetFormat(Format::EVALUATION);

    // acc = dct * ek (matrix product);
    // uses in-place * operators for the last call to dct[i] to gain performance improvement
    const std::vector<std::vector<NativePoly>>& ev = ek->GetElements();
    acc->GetElements()[0]                          = (dct[0] * ev[0][0]);
    for (uint32_t l = 1; l < digitsG2; ++l)
        acc->GetElements()[0] += (dct[l] * ev[l][0]);
    acc->GetElements()[1] = (dct[0] *= ev[0][1]);
    for (uint32_t l = 1; l < digitsG2; ++l)
        acc->GetElements()[1] += (dct[l] *= ev[l][1]);
}

// */

// TODO: remove this function here and in .h file
// Encryption as described in Section 5 of https://eprint.iacr.org/2014/816
// skNTT corresponds to the secret key z
RingGSWEvalKey RingGSWAccumulatorXZDDF::KeyGenDM(const std::shared_ptr<RingGSWCryptoParams>& params,
                                              const NativePoly& skNTT, LWEPlaintext m) const {
    const auto& Gpow       = params->GetGPower();
    const auto& polyParams = params->GetPolyParams();

    DiscreteUniformGeneratorImpl<NativeVector> dug;
    NativeInteger Q{params->GetQ()};
    dug.SetModulus(Q);

    // Reduce mod q (dealing with negative number as well)
    uint64_t q = params->Getq().ConvertToInt();
    uint32_t N = params->GetN();
    int64_t mm = (((m % q) + q) % q) * (2 * N / q);
    bool isReducedMM;
    if ((isReducedMM = (mm >= N)))
        mm -= N;

    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG2{(params->GetDigitsG() - 1) << 1};
    std::vector<NativePoly> tempA(digitsG2, NativePoly(dug, polyParams, Format::COEFFICIENT));
    RingGSWEvalKeyImpl result(digitsG2, 2);

    for (uint32_t i = 0; i < digitsG2; ++i) {
        result[i][0] = tempA[i];
        tempA[i].SetFormat(Format::EVALUATION);
        result[i][1] = NativePoly(params->GetDgg(), polyParams, Format::COEFFICIENT);
        if (!isReducedMM)
            result[i][i & 0x1][mm].ModAddFastEq(Gpow[(i >> 1) + 1], Q);
        else
            result[i][i & 0x1][mm].ModSubFastEq(Gpow[(i >> 1) + 1], Q);
        result[i][0].SetFormat(Format::EVALUATION);
        result[i][1].SetFormat(Format::EVALUATION);
        result[i][1] += (tempA[i] *= skNTT);
    }
    return std::make_shared<RingGSWEvalKeyImpl>(result);
}



/////////


// /*

// Key generation as described in Section 4 of https://eprint.iacr.org/2014/816
RingGSWACCKey RingGSWAccumulatorXZDDF::KeyGenAcc(const std::shared_ptr<RingGSWCryptoParams>& params,
                                              const NativePoly& skNTT, ConstLWEPrivateKey& LWEsk) const {
    auto sv = LWEsk->GetElement();
    auto q = params->Getq().ConvertToInt();
    auto N = params->GetN();
    if ((N % q) != 0)
        OPENFHE_THROW(config_error, "q does not divide N");        // which error type message to use?    
    if (!skNTT.InverseExists()) {
        OPENFHE_THROW(math_error, "GSW secret key (polynomial) has no inverse");        // which error type message to use?
    }
    
    uint32_t n(sv.GetLength());
    std::cout << "n: " << n << std::endl; 
    RingGSWACCKey ek = std::make_shared<RingGSWACCKeyImpl>(1, 1, (n+1) + (q-1));

    // evk_0
    auto s0{sv[0].ConvertToInt<int32_t>()};
    auto xExpS0 = GetXPower(params, s0);
    xExpS0.SetFormat(Format::EVALUATION);
    (*ek)[0][0][0] = NTRUencrypt(params, skNTT, xExpS0 * skNTT.MultiplicativeInverse());

/* 
    std::cout << "evk00 first: " << std::endl;
    auto p = (*ek)[0][0][0]->GetElements()[0][0];
    auto values = p.GetValues();
    auto len = values.GetLength();
    std::cout << "Length: " << len << std::endl;
    std::cout << "Elements:";
    for (long unsigned int i = 0; i < len; i++) {
        std::cout << " " << values[i].ConvertToInt();
    }
    std::cout << std::endl;
    int tmp; std::cin >> tmp;
 */

    // evk_i, i=1..(n-1)
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(n))
    for (uint32_t i = 1; i < n; ++i) {
        auto s{sv[i].ConvertToInt<int32_t>()};
        (*ek)[0][0][i] = NTRUencrypt(params, skNTT, GetXPower(params, s));
    }

    // evk_n
    auto negSum = -sv[0].ConvertToInt<int32_t>();
    for (uint32_t i = 1; i < n; ++i) {
        negSum -= sv[i].ConvertToInt<int32_t>();
    }
    (*ek)[0][0][n] = NTRUencrypt(params, skNTT, GetXPower(params, negSum));

    // ksk_j
    for (uint32_t i = 1; i < q; ++i) {
        auto autotrans = skNTT.AutomorphismTransform((2*N/q)*i + 1);
        autotrans.SetFormat(Format::EVALUATION);
        auto inv = skNTT.MultiplicativeInverse();
        (*ek)[0][0][n+i] = NTRUencrypt(params, skNTT, autotrans*inv);
    }

    return ek;
}
// */

// Vector NTRU encryption as described in Section 3 of https://eprint.iacr.org/2023/1564.pdf
// skNTT corresponds to the secret key f
RingGSWEvalKey RingGSWAccumulatorXZDDF::NTRUencrypt(const std::shared_ptr<RingGSWCryptoParams>& params,
                                        const NativePoly& skNTT, NativePoly m) const {

    const int64_t tau = 1; // TODO: add to RGSW params, and check if int64_t is best size

    const auto& powersG = params->GetGPower();
    std::vector<NativePoly> noisePols(powersG.size(), NativePoly(params->GetDgg(), params->GetPolyParams(), Format::COEFFICIENT));
    RingGSWEvalKeyImpl result(1, powersG.size());

    for (uint32_t i = 0; i < powersG.size(); ++i) {
        result[0][i] = noisePols[i].Times(tau);
        result[0][i].SetFormat(Format::EVALUATION);
        result[0][i] *= skNTT.MultiplicativeInverse();      // TODO: try sending inverse instead;  TODO: add check that skNTT is in evaluation mode??
        result[0][i].SetFormat(Format::COEFFICIENT);
        m.SetFormat(Format::COEFFICIENT);
        result[0][i] += m.Times(powersG[i]);             // TODO: ModAddFastEq instead? Need extra loop though
    }
    return std::make_shared<RingGSWEvalKeyImpl>(result);  // TODO: make_shared here or outside?
}

NativePoly RingGSWAccumulatorXZDDF::GetXPower(const std::shared_ptr<RingGSWCryptoParams>& params,
                                              LWEPlaintext m) const { // TODO: add format parameter with COEFFICIENT as default?
    const auto& polyParams = params->GetPolyParams();
    auto cycOrd = polyParams->GetCyclotomicOrder();
    uint32_t N = params->GetN();
    int64_t mm = (((m % cycOrd) + cycOrd) % cycOrd);
    int sign = 1;
    if (mm >= N) {
        mm -= N;
        sign = -1;
    }
    NativePoly mPol(polyParams, Format::COEFFICIENT, true);
    mPol[mm].SetValue(NativeInteger("1"));
    mPol = mPol.Times(sign);
    return mPol;
}

NativePoly RingGSWAccumulatorXZDDF::GetRotPol(const std::shared_ptr<RingGSWCryptoParams>& params, uint32_t expTransform = 1) const {  // TODO: is uint32_t good here?
    NativePoly r(params->GetPolyParams(), Format::COEFFICIENT, true);

    for (uint32_t i = 0; i < (params->Getq().ConvertToInt() >> 1); ++i) {
        r += GetXPower(params, -i*expTransform).Times(i);
    
/*         std::cout << i << std::endl;
        auto values = r.GetValues();
        auto len = values.GetLength();
        std::cout << "Length: " << len << std::endl;
        std::cout << "Elements:";
        for (long unsigned int i = 0; i < len; i++) {
            std::cout << " " << values[i].ConvertToInt();
        }
        std::cout << std::endl;
        int tmp; std::cin >> tmp; */
    }    
    return r;
}


NativePoly RingGSWAccumulatorXZDDF::ExternProd(const std::shared_ptr<RingGSWCryptoParams>& params, NativePoly c1, std::vector<PolyImpl<NativeVector>> c2) const {
    NativePoly extProd(params->GetPolyParams(), Format::COEFFICIENT, true);

    int baseBits = 0;
    uint32_t baseGtmp = params->GetBaseG();
    while (baseGtmp >>= 1) ++baseBits;              // TODO: move baseBits to params

    c1.SetFormat(Format::COEFFICIENT);
    auto bitdecom = c1.BaseDecompose(baseBits, true);

    for (unsigned int i = 0; i < bitdecom.size(); ++i) {
        auto p1 = bitdecom[i];
        auto p2 = c2[i];
        p1.SetFormat(Format::EVALUATION);
        p2.SetFormat(Format::EVALUATION);
        auto prod = p1 * p2;
        prod.SetFormat(Format::COEFFICIENT);
        extProd += prod;
    }

    return extProd;
}


// /*
void RingGSWAccumulatorXZDDF::EvalAcc(const std::shared_ptr<RingGSWCryptoParams>& params, ConstRingGSWACCKey& ek,
                                   RLWECiphertext& acc, const NativeVector& a) const {

    const int64_t delta = 1; // TODO: add to RGSW params, and check if int64_t is best size

    auto N = params->GetN();
    auto q       = params->Getq().ConvertToInt();
    uint32_t n = a.GetLength();                     // TODO: add -1 when b is appended to end of a
    std::cout << "N: " << N << ", q: " << q << ", n: " << n << std::endl;
    std::cout << "B: " << params->GetBaseG() << std::endl;

    NativeInteger wNow = 2*(N/q)*a[0].ConvertToInt() + 1;
    NativeInteger wNext;
    std::cout << "w0: " << wNow.ConvertToInt() << std::endl;

    auto b = acc->GetElements()[1][0].ConvertToInt();
    std::cout << "b2: " << b << std::endl;

    auto tmpRotPolExp = 2*(N/q)*wNow.ModInverse(2*N).ConvertToInt();
    std::cout << "rot pol exp: " << tmpRotPolExp << std::endl;


    auto accPol = GetRotPol(params, 2*(N/q)*wNow.ModInverse(2*N).ConvertToInt());
    auto xPower0 = GetXPower(params, -2*(N/q)*b*wNow.ModInverse(2*N).ConvertToInt());
    accPol.SetFormat(Format::EVALUATION);
    xPower0.SetFormat(Format::EVALUATION);
    accPol = accPol * xPower0;
    accPol.SetFormat(Format::COEFFICIENT);
    accPol = accPol.Times(delta);   

    for (uint32_t i = 0; i < n; ++i) {
        auto ev = (*ek)[0][0][i]->GetElements()[0];
        // std::cout << ev.GetValues()[0].ConvertToInt() << " ";
        
        accPol = ExternProd(params, accPol, ev);

        if (i < a.GetLength() - 1) {
            wNext = 2*(N/q)*a[i+1].ConvertToInt() + 1;
        }
        else {
            wNext = 1;
        }
        
        auto wProd = (wNow * wNext.ModInverse(2*N)).ConvertToInt() % (2*N);
        if (wProd != 1) {
            int kskInd = (wProd-1)*q/(2*N) + n;
            auto ksk = (*ek)[0][0][kskInd]->GetElements()[0];
            accPol = ExternProd(params, accPol.AutomorphismTransform(wProd), ksk);
        }

    }
    auto ev = (*ek)[0][0][n]->GetElements()[0];
    accPol = ExternProd(params, accPol, ev);

    std::cout << "evalAcc done" << std::endl;
    acc->GetElements()[0] = accPol;

}

// AP Accumulation as described in https://eprint.iacr.org/2020/086
void RingGSWAccumulatorXZDDF::AddToAccDM(const std::shared_ptr<RingGSWCryptoParams>& params, ConstRingGSWEvalKey& ek,
                                      RLWECiphertext& acc) const {
    std::vector<NativePoly> ct(acc->GetElements());
    ct[0].SetFormat(Format::COEFFICIENT);
    ct[1].SetFormat(Format::COEFFICIENT);

    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG2{(params->GetDigitsG() - 1) << 1};
    std::vector<NativePoly> dct(digitsG2, NativePoly(params->GetPolyParams(), Format::COEFFICIENT, true));

    SignedDigitDecompose(params, ct, dct);

#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(digitsG2))
    for (uint32_t j = 0; j < digitsG2; ++j)
        dct[j].SetFormat(Format::EVALUATION);

    // acc = dct * ek (matrix product);
    // uses in-place * operators for the last call to dct[i] to gain performance improvement
    const std::vector<std::vector<NativePoly>>& ev = ek->GetElements();
    acc->GetElements()[0]                          = (dct[0] * ev[0][0]);
    for (uint32_t l = 1; l < digitsG2; ++l)
        acc->GetElements()[0] += (dct[l] * ev[l][0]);
    acc->GetElements()[1] = (dct[0] *= ev[0][1]);
    for (uint32_t l = 1; l < digitsG2; ++l)
        acc->GetElements()[1] += (dct[l] *= ev[l][1]);
}
// */
};  // namespace lbcrypto
