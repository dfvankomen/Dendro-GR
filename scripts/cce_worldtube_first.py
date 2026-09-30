##
# @brief: FIRST-ORDER shift-driver worldtube config: supplies Dendro-GR's
#          natively-evolved auxiliary shift B^i directly to SpECTRE CCE's
#          `PreprocessCceWorldtube` (AdmMetricNodal, first-order driver),
#          and emits a ready-to-use PreprocessCceWorldtube.yaml.
# @date: 2026-09-16
#
# Context: SpECTRE's `PreprocessCceWorldtube` (AdmMetricNodal) supports two
# mutually-exclusive ways to reconstruct the shift's time derivative -- see
# cce_worldtube_second.py for the sibling second-order (Gamma-driver eta)
# path. THIS script covers the first-order path, which reconstructs
#     d(beta^i)/dt = FirstOrderDriverFactor * B^i + (advection)
# reading B^i directly from the `/AuxiliaryShift` dataset in the worldtube
# HDF5 file -- SpECTRE does NOT compute B^i itself for this path, unlike
# the eta-based reconstruction in the second-order path.
#
# Dendro-GR already evolves B^i natively as one of its 24 BSSN fields
# (U_B0, U_B1, U_B2, see BSSN_GR/include/grDef.h) -- it is not derived or
# assumed, just interpolated onto the worldtube sphere like any other
# field (via the same ot::da::interpolateToCoords call used everywhere
# else in Dendro-GR, e.g. the BHaHAHA apparent-horizon integration). This
# script therefore does NOT need to evaluate any analytic prescription --
# unlike cce_worldtube_second.py's eta(r) evaluation, there is nothing to
# compute here beyond confirming/echoing the correct FirstOrderDriverFactor
# constant and emitting the YAML. The actual /AuxiliaryShift dataset must
# be populated separately, by the (currently on-hold, see Dendro_CCE
# writeup Section 7 step 4) worldtube-writer module that interpolates
# Dendro-GR's evolved B^i field onto the worldtube's angular grid.
#
# FirstOrderDriverFactor and Dendro-GR's own compiled-in prefactor:
# Dendro-GR's shift RHS (BSSN_GR/src/bssneqs_eta_const_standard_gauge.cpp)
# is
#     d(beta^i)/dt = DENDRO_1 * B^i + (advection),
#     DENDRO_1 = (3/4) * alpha * BSSN_LAMBDA_F[1] + (3/4) * BSSN_LAMBDA_F[0]
# With the compiled-in default BSSN_LAMBDA_F = [1.0, 0.0] (BSSN_GR/src/
# parameters.cpp, and NOT overridden in the sample q1/q4/q8 .par.toml
# configs), DENDRO_1 = 3/4 exactly -- a true constant, independent of
# alpha -- which is exactly SpECTRE's own example FirstOrderDriverFactor
# value ("0.75 # Eq. 4.89 B&S", Baumgarte & Shapiro, "Numerical Relativity:
# Solving Einstein's Equations on the Computer," Cambridge Univ. Press
# 2010, Eq. 4.89 -- see spectre/docs/References.bib entry `BaumgarteShapiro`
# and spectre/docs/Tutorials/CCE.md).
#
# CAVEAT (kept honest, by analogy with the eta script): if BSSN_LAMBDA_F[1]
# is nonzero (the alpha-weighted "shock-avoiding" Gamma-driver variant),
# DENDRO_1 becomes alpha-dependent -- i.e. it varies pointwise over the
# worldtube sphere, and can no longer be represented as the single YAML
# constant FirstOrderDriverFactor expects. This script warns loudly and
# refuses to silently emit a wrong constant in that case; the sample
# q1/q4/q8 configs all use the default (BSSN_LAMBDA_F[1] = 0), where this
# is a non-issue.
#
# NOTE: whether the /ConformalChristoffel (Gt^i) side of the *second-order*
# path can be made to work is still under review by a group member -- this
# first-order script is the currently-recommended default in the meantime
# (see Dendro_CCE_v1.3.md, Section 6.2), since it requires no analytic
# eta-matching and no invented ConformalChristoffelFactor constant, only
# Dendro-GR's own natively-evolved B^i field.

import argparse
import sys

# compiled-in default, BSSN_GR/src/parameters.cpp: BSSN_LAMBDA_F = {1.0, 0.0}
DEFAULT_LAMBDA_F0 = 1.0
DEFAULT_LAMBDA_F1 = 0.0
# Dendro-GR's compiled-in prefactor with the defaults above (see docstring)
DEFAULT_FIRST_ORDER_DRIVER_FACTOR = 0.75


def first_order_driver_factor(lambda_f0=DEFAULT_LAMBDA_F0,
                               lambda_f1=DEFAULT_LAMBDA_F1,
                               alpha=None):
    '''
    @brief Computes Dendro-GR's compiled-in shift-driver prefactor
           DENDRO_1 = (3/4)*alpha*lambda_f1 + (3/4)*lambda_f0, exactly as it
           appears in BSSN_GR/src/bssneqs_eta_const_standard_gauge.cpp.
           Raises ValueError if lambda_f1 != 0 and no representative alpha
           is given, since the factor is then pointwise-varying and cannot
           be represented as a single constant -- see module docstring.
    '''
    if lambda_f1 != 0.0:
        if alpha is None:
            raise ValueError(
                'BSSN_LAMBDA_F[1] = %g (nonzero, the alpha-weighted '
                '"shock-avoiding" Gamma-driver variant) -- the shift-driver '
                'prefactor is then alpha-dependent (varies pointwise over '
                'the worldtube sphere) and CANNOT be represented exactly as '
                'a single FirstOrderDriverFactor constant. Pass --alpha '
                'with a representative lapse value if you want an '
                'approximate constant anyway (not recommended -- consider '
                'the second-order path, cce_worldtube_second.py, or a '
                'worldtube-writer-side per-point implementation instead).'
                % lambda_f1)
        print('WARNING: BSSN_LAMBDA_F[1] = %g is nonzero -- the shift-driver '
              'prefactor is alpha-dependent, not a true constant. Using '
              'alpha=%g as a representative value; this is an '
              'APPROXIMATION, not exact, unlike the lambda_f1=0 case.'
              % (lambda_f1, alpha), file=sys.stderr)
    else:
        alpha = 0.0  # unused when lambda_f1 == 0
    return 0.75 * alpha * lambda_f1 + 0.75 * lambda_f0


YAML_TEMPLATE = """\
# Generated by cce_worldtube_first.py -- Dendro-GR -> SpECTRE CCE worldtube
# preprocessing config, FIRST-ORDER shift driver. Mirrors SpECTRE's own
# tests/InputFiles/PreprocessCceWorldtube/AdmFirstOrderDriverPreprocessCceWorldtube.yaml.
#
# IMPORTANT: this format requires a `/AuxiliaryShift` dataset (Dendro-GR's
# natively-evolved B^i field, U_B0/U_B1/U_B2, interpolated onto the
# worldtube's angular grid) in the input worldtube HDF5 file -- unlike the
# second-order path, SpECTRE does not compute B^i itself here.

InputH5File: {input_h5}
OutputH5File: {output_h5}
InputDataFormat:
  AdmMetricNodal:
    Lapse:
      Advective: True
    Shift:
      Advective: True
      FirstOrderDriverFactor: {driver_factor:.10g} # Eq. 4.89 B&S (Baumgarte & Shapiro 2010); Dendro-GR compiled-in value with BSSN_LAMBDA_F=[{lambda_f0:g},{lambda_f1:g}]
ExtractionRadius: {radius:g}
FixSpecNormalization: False
DescendingM: False
BufferDepth: Auto
LMaxFactor: {lmax_factor:g}
"""


def main():
    parser = argparse.ArgumentParser(
        description='Emit a ready-to-use PreprocessCceWorldtube.yaml for '
                     'the first-order shift-driver path, using Dendro-GR\'s '
                     'compiled-in shift-RHS prefactor as FirstOrderDriverFactor. '
                     'Unlike cce_worldtube_second.py, this does not need to '
                     'evaluate any per-radius analytic formula -- B^i is '
                     'supplied directly from Dendro-GR\'s own evolved field.')
    parser.add_argument('--radius', type=float, required=True,
                         help='Worldtube/extraction coordinate radius, in '
                              'units of total mass M.')
    parser.add_argument('--lambda-f0', type=float, default=DEFAULT_LAMBDA_F0,
                         help='Override BSSN_LAMBDA_F[0] (default: compiled-in '
                              'Dendro-GR default, %g).' % DEFAULT_LAMBDA_F0)
    parser.add_argument('--lambda-f1', type=float, default=DEFAULT_LAMBDA_F1,
                         help='Override BSSN_LAMBDA_F[1] (default: compiled-in '
                              'Dendro-GR default, %g). If nonzero, --alpha is '
                              'required (see module docstring caveat).'
                              % DEFAULT_LAMBDA_F1)
    parser.add_argument('--alpha', type=float, default=None,
                         help='Representative lapse value, only used (and '
                              'only an approximation) if --lambda-f1 is '
                              'nonzero.')
    parser.add_argument('--write-yaml', type=str, default=None,
                         help='If given, write a full PreprocessCceWorldtube.yaml '
                              'to this path with the driver factor plugged in.')
    parser.add_argument('--input-h5', type=str, default='InputFilename.h5',
                         help='InputH5File value for --write-yaml (default: '
                              'placeholder, edit before use).')
    parser.add_argument('--output-h5', type=str, default=None,
                         help='OutputH5File value for --write-yaml (default: '
                              'derived from --radius, e.g. ReducedWorldtubeR0100.h5).')
    parser.add_argument('--lmax-factor', type=float, default=3,
                         help='LMaxFactor value for --write-yaml (default: 3).')
    args = parser.parse_args()

    try:
        driver_factor = first_order_driver_factor(
            args.lambda_f0, args.lambda_f1, args.alpha)
    except ValueError as exc:
        print('ERROR: %s' % exc, file=sys.stderr)
        return 1

    print('BSSN_LAMBDA_F = [%g, %g]' % (args.lambda_f0, args.lambda_f1))
    print('FirstOrderDriverFactor = %.10g' % driver_factor)

    if args.write_yaml:
        output_h5 = args.output_h5 or (
            'ReducedWorldtubeR%04d.h5' % round(args.radius))
        yaml_text = YAML_TEMPLATE.format(
            input_h5=args.input_h5,
            output_h5=output_h5,
            driver_factor=driver_factor,
            lambda_f0=args.lambda_f0,
            lambda_f1=args.lambda_f1,
            radius=args.radius,
            lmax_factor=args.lmax_factor,
        )
        with open(args.write_yaml, 'w') as f:
            f.write(yaml_text)
        print('Wrote %s' % args.write_yaml)
        print('REMINDER: populate /AuxiliaryShift in %s with Dendro-GR\'s '
              'interpolated B^i field before running PreprocessCceWorldtube.'
              % args.input_h5)

    return 0


if __name__ == '__main__':
    sys.exit(main())
