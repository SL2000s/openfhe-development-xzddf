// TODO: now copy of AP/DM, change to XZDDF

#include "rgsw-acc-xzddf.h"

#include <string>

namespace lbcrypto {

// Key generation as described in Section 4 of https://eprint.iacr.org/2023/1564
RingGSWACCKey RingGSWAccumulatorXZDDF::KeyGenAcc(const std::shared_ptr<RingGSWCryptoParams>& params,
                                              const NativePoly& skNTT, ConstLWEPrivateKey& LWEsk) const {
    auto sv = LWEsk->GetElement();
    auto q = params->Getq().ConvertToInt();
    auto N = params->GetN();
    if ((N % q) != 0)
        OPENFHE_THROW(config_error, "q does not divide N");                             // TODO: which error type message to use?    
    if (!skNTT.InverseExists()) {
        OPENFHE_THROW(math_error, "GSW secret key (polynomial) has no inverse");        // TODO: which error type message to use?
    }
    uint32_t n = sv.GetLength();
    std::cout << "n: " << n << std::endl; 
    
    RingGSWACCKey ek = std::make_shared<RingGSWACCKeyImpl>(1, 1, (n+1) + (q-1));
    const auto skNTTinv = skNTT.MultiplicativeInverse();                            // will be on evaluation form

    // evk_0
    auto s0 = sv[0].ConvertToInt<int32_t>();
    auto xExpS0 = GetXPower(params, s0);
    xExpS0.SetFormat(Format::EVALUATION);
    // skNTTinv.SetFormat(Format::EVALUATION);
    (*ek)[0][0][0] = NTRUencrypt(params, skNTTinv, xExpS0 * skNTTinv);

    // evk_i, i=1..(n-1)
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(n))
    for (uint32_t i = 1; i < n; ++i) {
        auto s = sv[i].ConvertToInt<int32_t>();
        (*ek)[0][0][i] = NTRUencrypt(params, skNTTinv, GetXPower(params, s));
    }

    // evk_n
    auto negSum = -sv[0].ConvertToInt<int32_t>();
    for (uint32_t i = 1; i < n; ++i) {
        negSum -= sv[i].ConvertToInt<int32_t>();
    }
    (*ek)[0][0][n] = NTRUencrypt(params, skNTTinv, GetXPower(params, negSum));

    // ksk_j
    for (uint32_t i = 1; i < q; ++i) {
        auto autotrans = skNTT.AutomorphismTransform((2*N/q)*i + 1);
        autotrans.SetFormat(Format::EVALUATION);
        // skNTTinv.SetFormat(Format::EVALUATION);
        (*ek)[0][0][n+i] = NTRUencrypt(params, skNTTinv, autotrans*skNTTinv);
    }

    return ek;
}

// Vector NTRU encryption as described in Section 3 of https://eprint.iacr.org/2023/1564.pdf
// skNTT corresponds to the secret key f
RingGSWEvalKey RingGSWAccumulatorXZDDF::NTRUencrypt(const std::shared_ptr<RingGSWCryptoParams>& params,
                                        const NativePoly& skNTTinv, NativePoly m) const {

    const int64_t tau = 1;                                  // TODO: add to RGSW params, and check if int64_t is best size

    const auto& powersG = params->GetGPower();
    auto d = params->GetDigitsG();
    std::vector<NativePoly> noisePols(d, NativePoly(params->GetDgg(), params->GetPolyParams(), Format::COEFFICIENT));
    RingGSWEvalKeyImpl result(1, d);
    for (uint32_t i = 0; i < d; ++i) {
        result[0][i] = noisePols[i].Times(tau);
        result[0][i].SetFormat(Format::EVALUATION);
        result[0][i] *= skNTTinv;                           // TODO: add check that skNTT is in evaluation mode??
        result[0][i].SetFormat(Format::COEFFICIENT);
        m.SetFormat(Format::COEFFICIENT);
        result[0][i] += m.Times(powersG[i]);                // TODO: ModAddFastEq instead? Need extra loop though
    }
    return std::make_shared<RingGSWEvalKeyImpl>(result);    // TODO: make_shared here or outside?
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

NativePoly RingGSWAccumulatorXZDDF::GetRotPol(const std::shared_ptr<RingGSWCryptoParams>& params, uint32_t expTransform = 1) const {  // TODO: is uint32_t good here?  TODO: =1 -> .h file
    NativePoly r(params->GetPolyParams(), Format::COEFFICIENT, true);
    auto q = params->Getq().ConvertToInt();
    for (uint32_t i = 0; i < (q >> 1); ++i) {
        // r += GetXPower(params, -i*expTransform).Times( ((i+ (q>>3)) / (q>>2)) % 4 ); // r(x) = sum_{i=1..(q-1)}[round(i/(q/t))*X^(-i)]
        
        int64_t coeff = 1;
        if (i >= (q>>2)) {
            coeff = -1;
        }
        r += GetXPower(params, -i*expTransform).Times(coeff);       // TODO: optimize this -- times(1) unnecessary



        // std::cout << " rot pol:   i: " << i << " q: " << q << " prod: " << ((i+ (q>>3)) / (q>>2)) % 4  << std::endl;
        // r += GetXPower(params, -i*expTransform).Times(i);
    
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

/*     auto values = r.GetValues();
    auto len = values.GetLength();
    std::cout << "Length: " << len << std::endl;
    std::cout << "Elements:";
    for (long unsigned int i = 0; i < len; i++) {
        std::cout << " " << values[i].ConvertToInt();
    }
    std::cout << std::endl;
    // int tmp; std::cin >> tmp;
 */

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

void RingGSWAccumulatorXZDDF::EvalAcc(const std::shared_ptr<RingGSWCryptoParams>& params, ConstRingGSWACCKey& ek,
                                   RLWECiphertext& acc, const NativeVector& a) const {

    const int64_t delta = 7*(params->GetQ().ConvertToInt() >> 3); // TODO: add to RGSW params, and check if int64_t is best size
    // const int64_t delta = params->GetQ().ConvertToInt() >> 2; // TODO: add to RGSW params, and check if int64_t is best size

    auto N = params->GetN();
    auto q       = params->Getq().ConvertToInt();
    uint32_t n = a.GetLength();                     // TODO: add -1 when b is appended to end of a
    std::cout << "N: " << N << ", q: " << q << ", n: " << n << " Q: " << params->GetQ().ConvertToInt() << std::endl;
    std::cout << "B: " << params->GetBaseG() << std::endl;

    NativeInteger wNow = ( (2*2*N - 2*(N/q)*a[0].ConvertToInt()) + 1) % (2*N);  // w_i = (-(2*N/q)*a_i+1)

    NativeInteger wNext;
    std::cout << "w0: " << wNow.ConvertToInt() << std::endl;

    auto b = acc->GetElements()[1][0].ConvertToInt();
    std::cout << "b2: " << b << std::endl;

    auto tmpRotPolExp = 2*(N/q)*wNow.ModInverse(2*N).ConvertToInt();
    std::cout << "rot pol exp: " << tmpRotPolExp << std::endl;

    auto accPol = GetRotPol(params, 2*(N/q)*wNow.ModInverse(2*N).ConvertToInt());
    auto xPower0 = GetXPower(params, 2*(N/q)*(b)*wNow.ModInverse(2*N).ConvertToInt());   // q>>3: to transform 
    // auto xPower0 = GetXPower(params, 2*(N/q)*b*wNow.ModInverse(2*N).ConvertToInt());
    accPol.SetFormat(Format::EVALUATION);
    xPower0.SetFormat(Format::EVALUATION);
    accPol = accPol * xPower0;
    accPol.SetFormat(Format::COEFFICIENT);
    accPol = accPol.Times(delta);   

    for (uint32_t i = 0; i < n; ++i) {
        auto ev = (*ek)[0][0][i]->GetElements()[0];
        
        accPol = ExternProd(params, accPol, ev);

        if (i < a.GetLength() - 1) {
            wNext = ((2*2*N - 2*(N/q)*a[i+1].ConvertToInt()) + 1) % (2*N);
        }
        else {
            wNext = 1;
        }
        
        auto wProd = (wNow * wNext.ModInverse(2*N)).ConvertToInt() % (2*N);
        if (wProd != 1) {
            int kskInd = (wProd-1)*q/(2*N) + n;
            auto ksk = (*ek)[0][0][kskInd]->GetElements()[0];
            accPol.SetFormat(Format::COEFFICIENT);                                      // TODO: check that automorphisms work for coefficient format
            accPol = ExternProd(params, accPol.AutomorphismTransform(wProd), ksk);
        }

        wNow = wNext;
    }
    auto ev = (*ek)[0][0][n]->GetElements()[0];
    accPol = ExternProd(params, accPol, ev);

    acc->GetElements()[0] = accPol;
}

};  // namespace lbcrypto
