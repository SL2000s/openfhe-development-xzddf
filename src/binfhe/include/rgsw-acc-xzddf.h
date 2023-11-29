// TODO: now copy of AP/DM, change to XZDDF

#ifndef _RGSW_ACC_XZDDF_H_
#define _RGSW_ACC_ZDDF_H_

#include "rgsw-acc.h"

#include <memory>

namespace lbcrypto {

/**
 * @brief Ring GSW accumulator schemes described in
 * https://eprint.iacr.org/2023/1564
 */
class RingGSWAccumulatorXZDDF final : public RingGSWAccumulator {
public:
    RingGSWAccumulatorXZDDF() = default;

    /**
   * Key generation for internal Ring GSW as described in https://eprint.iacr.org/2020/086
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param skNTT secret key polynomial in the EVALUATION representation
   * @param LWEsk the secret key
   * @return a shared pointer to the resulting keys
   */
    RingGSWACCKey KeyGenAcc(const std::shared_ptr<RingGSWCryptoParams>& params, const NativePoly& skNTT,
                            ConstLWEPrivateKey& LWEsk) const override;

    /**
   * Main accumulator function used in bootstrapping - AP variant
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param ek the accumulator key
   * @param acc previous value of the accumulator
   * @param a value to update the accumulator with
   */
    void EvalAcc(const std::shared_ptr<RingGSWCryptoParams>& params, ConstRingGSWACCKey& ek, RLWECiphertext& acc,
                 const NativeVector& a) const override;

private:
    // TODO: change explainations
    /**
   * DM Key generation for internal Ring GSW as described in https://eprint.iacr.org/2014/816
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param skNTT secret key polynomial in the EVALUATION representation
   * @param m a plaintext
   * @return a shared pointer to the resulting keys
   */
    RingGSWEvalKey NTRUencrypt(const std::shared_ptr<RingGSWCryptoParams>& params, const NativePoly& skNTTinv,
                                   NativePoly m) const;

    // TODO: add explainations parameters
    NativePoly GetXPower(const std::shared_ptr<RingGSWCryptoParams>& params, LWEPlaintext m) const;
    NativePoly GetRotPol(const std::shared_ptr<RingGSWCryptoParams>& params, uint32_t expTransform) const;
    NativePoly ExternProd(const std::shared_ptr<RingGSWCryptoParams> &params, NativePoly c1, std::vector<PolyImpl<NativeVector>> c2) const;
};

}  // namespace lbcrypto

#endif  // _RGSW_ACC_XZDDF_H_
