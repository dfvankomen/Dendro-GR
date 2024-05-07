// Codgen: generating unstage version
// Codgen: using standard gauge
// Codgen: using eta const damping
//  Dendro: {{{
//  Dendro: original ops: 623312
//  Dendro: printing temp variables
const double DENDRO_0 =
    (3.0 / 4.0) * alpha[pp] * lambda_f[1] + (3.0 / 4.0) * lambda_f[0];
const double DENDRO_1  = 2 * alpha[pp];
const double DENDRO_2  = (4.0 / 3.0) * grad_0_beta0[pp];
const double DENDRO_3  = 2 * gt1[pp];
const double DENDRO_4  = 2 * gt2[pp];
const double DENDRO_5  = (2.0 / 3.0) * gt0[pp];
const double DENDRO_6  = (1.0 / 3.0) * gt1[pp];
const double DENDRO_7  = (2.0 / 3.0) * grad_2_beta2[pp];
const double DENDRO_8  = (1.0 / 3.0) * gt2[pp];
const double DENDRO_9  = (2.0 / 3.0) * grad_1_beta1[pp];
const double DENDRO_10 = (2.0 / 3.0) * grad_0_beta0[pp];
const double DENDRO_11 = (4.0 / 3.0) * grad_1_beta1[pp];
const double DENDRO_12 = 2 * gt4[pp];
const double DENDRO_13 = (1.0 / 3.0) * gt4[pp];
const double DENDRO_14 = (4.0 / 3.0) * grad_2_beta2[pp];
const double DENDRO_15 = (2.0 / 3.0) * chi[pp];
const double DENDRO_16 = grad_0_beta0[pp] + grad_1_beta1[pp] + grad_2_beta2[pp];
const double DENDRO_17 = 2 * At1[pp];
const double DENDRO_18 = 2 * At2[pp];
const double DENDRO_19 = gt2[pp] * gt4[pp];
const double DENDRO_20 = -DENDRO_19 + gt1[pp] * gt5[pp];
const double DENDRO_21 = At0[pp] * DENDRO_20;
const double DENDRO_22 = gt0[pp] * gt4[pp] - gt1[pp] * gt2[pp];
const double DENDRO_23 = At2[pp] * DENDRO_22;
const double DENDRO_24 = pow(gt2[pp], 2);
const double DENDRO_25 = -DENDRO_24 + gt0[pp] * gt5[pp];
const double DENDRO_26 = pow(gt4[pp], 2);
const double DENDRO_27 = pow(gt1[pp], 2);
const double DENDRO_28 = gt3[pp] * gt5[pp];
const double DENDRO_29 = -DENDRO_19 * DENDRO_3 + DENDRO_24 * gt3[pp] +
                         DENDRO_26 * gt0[pp] + DENDRO_27 * gt5[pp] -
                         DENDRO_28 * gt0[pp];
const double DENDRO_30 = 1.0 / DENDRO_29;
const double DENDRO_31 = DENDRO_17 * DENDRO_30;
const double DENDRO_32 = gt1[pp] * gt4[pp] - gt2[pp] * gt3[pp];
const double DENDRO_33 = At2[pp] * DENDRO_32;
const double DENDRO_34 = -DENDRO_26 + DENDRO_28;
const double DENDRO_35 = At0[pp] * DENDRO_34;
const double DENDRO_36 = 2 * DENDRO_30;
const double DENDRO_37 = At0[pp] * DENDRO_36;
const double DENDRO_38 = -DENDRO_27 + gt0[pp] * gt3[pp];
const double DENDRO_39 = At2[pp] * DENDRO_38;
const double DENDRO_40 = DENDRO_18 * DENDRO_30;
const double DENDRO_41 =
    DENDRO_20 * grad_0_chi[pp] + DENDRO_22 * grad_2_chi[pp];
const double DENDRO_42 = -DENDRO_25 * grad_1_chi[pp] + DENDRO_41;
const double DENDRO_43 = 0.5 * DENDRO_42;
const double DENDRO_44 = 1.0 / chi[pp];
const double DENDRO_45 = DENDRO_44 * gt0[pp];
const double DENDRO_46 = 1.0 * grad_0_gt1[pp];
const double DENDRO_47 = 0.5 * grad_1_gt0[pp];
const double DENDRO_48 = DENDRO_46 - DENDRO_47;
const double DENDRO_49 = 1.0 * grad_0_gt2[pp];
const double DENDRO_50 = 0.5 * grad_2_gt0[pp];
const double DENDRO_51 = DENDRO_49 - DENDRO_50;
const double DENDRO_52 = 0.5 * grad_0_gt0[pp];
const double DENDRO_53 = DENDRO_20 * DENDRO_52 + DENDRO_22 * DENDRO_51;
const double DENDRO_54 = -DENDRO_25 * DENDRO_48 + DENDRO_53;
const double DENDRO_55 = DENDRO_43 * DENDRO_45 + DENDRO_54;
const double DENDRO_56 = 12 * grad_1_alpha[pp];
const double DENDRO_57 = DENDRO_30 * DENDRO_56;
const double DENDRO_58 = -DENDRO_22 * grad_1_chi[pp] +
                         DENDRO_32 * grad_0_chi[pp] +
                         DENDRO_38 * grad_2_chi[pp];
const double DENDRO_59 = -DENDRO_58;
const double DENDRO_60 =
    -DENDRO_22 * DENDRO_48 + DENDRO_32 * DENDRO_52 + DENDRO_38 * DENDRO_51;
const double DENDRO_61 = -0.5 * DENDRO_44 * DENDRO_59 * gt0[pp] + DENDRO_60;
const double DENDRO_62 = 12 * grad_2_alpha[pp];
const double DENDRO_63 = DENDRO_30 * DENDRO_62;
const double DENDRO_64 = DENDRO_34 * DENDRO_52;
const double DENDRO_65 = DENDRO_32 * DENDRO_51;
const double DENDRO_66 = DENDRO_20 * DENDRO_48;
const double DENDRO_67 = DENDRO_30 * DENDRO_66;
const double DENDRO_68 = 1.0 * grad_0_chi[pp];
const double DENDRO_69 = DENDRO_30 * gt0[pp];
const double DENDRO_70 = -DENDRO_20 * grad_1_chi[pp] +
                         DENDRO_32 * grad_2_chi[pp] +
                         DENDRO_34 * grad_0_chi[pp];
const double DENDRO_71 = -DENDRO_70;
const double DENDRO_72 = 0.5 * DENDRO_71;
const double DENDRO_73 = DENDRO_30 * DENDRO_64 + DENDRO_30 * DENDRO_65 +
                         DENDRO_44 * (DENDRO_68 - DENDRO_69 * DENDRO_72) -
                         DENDRO_67;
const double DENDRO_74  = 12 * grad_0_alpha[pp];
const double DENDRO_75  = pow(chi[pp], -2);
const double DENDRO_76  = pow(grad_0_chi[pp], 2);
const double DENDRO_77  = DENDRO_75 * DENDRO_76;
const double DENDRO_78  = 4.0 * DENDRO_30;
const double DENDRO_79  = DENDRO_20 * DENDRO_78;
const double DENDRO_80  = DENDRO_79 * grad2_0_1_gt0[pp];
const double DENDRO_81  = DENDRO_22 * DENDRO_78;
const double DENDRO_82  = DENDRO_81 * grad2_1_2_gt0[pp];
const double DENDRO_83  = pow(DENDRO_29, -2);
const double DENDRO_84  = DENDRO_20 * grad_0_gt3[pp];
const double DENDRO_85  = DENDRO_34 * grad_1_gt0[pp];
const double DENDRO_86  = grad_0_gt4[pp] + grad_1_gt2[pp] - grad_2_gt1[pp];
const double DENDRO_87  = DENDRO_32 * DENDRO_86;
const double DENDRO_88  = -DENDRO_84 + DENDRO_85 + DENDRO_87;
const double DENDRO_89  = DENDRO_32 * grad_0_gt5[pp];
const double DENDRO_90  = DENDRO_34 * grad_2_gt0[pp];
const double DENDRO_91  = grad_0_gt4[pp] - grad_1_gt2[pp] + grad_2_gt1[pp];
const double DENDRO_92  = DENDRO_20 * DENDRO_91;
const double DENDRO_93  = DENDRO_89 + DENDRO_90 - DENDRO_92;
const double DENDRO_94  = -DENDRO_20 * DENDRO_48 + DENDRO_64 + DENDRO_65;
const double DENDRO_95  = DENDRO_25 * grad_0_gt3[pp];
const double DENDRO_96  = DENDRO_20 * grad_1_gt0[pp];
const double DENDRO_97  = DENDRO_22 * DENDRO_86;
const double DENDRO_98  = DENDRO_96 + DENDRO_97;
const double DENDRO_99  = -DENDRO_95 + DENDRO_98;
const double DENDRO_100 = 0.25 * grad_0_gt3[pp];
const double DENDRO_101 = 1.0 * grad_1_gt1[pp];
const double DENDRO_102 = -DENDRO_101;
const double DENDRO_103 = DENDRO_25 * DENDRO_83;
const double DENDRO_104 = 4 * DENDRO_103;
const double DENDRO_105 = DENDRO_104 * DENDRO_99 * (-DENDRO_100 - DENDRO_102);
const double DENDRO_106 = DENDRO_32 * grad_2_gt0[pp];
const double DENDRO_107 = DENDRO_38 * grad_0_gt5[pp];
const double DENDRO_108 = DENDRO_106 + DENDRO_107 - DENDRO_22 * DENDRO_91;
const double DENDRO_109 = 0.25 * grad_0_gt5[pp];
const double DENDRO_110 = 1.0 * grad_2_gt2[pp];
const double DENDRO_111 = -DENDRO_110;
const double DENDRO_112 = DENDRO_109 + DENDRO_111;
const double DENDRO_113 = DENDRO_38 * DENDRO_83;
const double DENDRO_114 = 4 * DENDRO_113;
const double DENDRO_115 =
    DENDRO_20 * grad_2_gt0[pp] + DENDRO_22 * grad_0_gt5[pp];
const double DENDRO_116 = DENDRO_115 - DENDRO_25 * DENDRO_91;
const double DENDRO_117 = 0.25 * grad_0_gt4[pp];
const double DENDRO_118 = -DENDRO_117;
const double DENDRO_119 = 0.25 * grad_1_gt2[pp];
const double DENDRO_120 = 0.75 * grad_2_gt1[pp];
const double DENDRO_121 =
    DENDRO_114 * DENDRO_116 * (DENDRO_118 + DENDRO_119 + DENDRO_120);
const double DENDRO_122 = 0.75 * grad_1_gt2[pp];
const double DENDRO_123 = 0.25 * grad_2_gt1[pp];
const double DENDRO_124 = DENDRO_118 + DENDRO_122 + DENDRO_123;
const double DENDRO_125 = DENDRO_32 * grad_1_gt0[pp];
const double DENDRO_126 = DENDRO_38 * DENDRO_86;
const double DENDRO_127 = DENDRO_125 + DENDRO_126 - DENDRO_22 * grad_0_gt3[pp];
const double DENDRO_128 = DENDRO_34 * DENDRO_54;
const double DENDRO_129 = 4 * DENDRO_83;
const double DENDRO_130 = DENDRO_128 * DENDRO_129 * (DENDRO_46 + DENDRO_47);
const double DENDRO_131 = DENDRO_49 + DENDRO_50;
const double DENDRO_132 = 0.25 * grad_1_gt0[pp];
const double DENDRO_133 = DENDRO_132 * DENDRO_93;
const double DENDRO_134 = DENDRO_22 * DENDRO_83;
const double DENDRO_135 = 4 * DENDRO_134;
const double DENDRO_136 = 0.25 * grad_2_gt0[pp];
const double DENDRO_137 = DENDRO_136 * DENDRO_88;
const double DENDRO_138 = DENDRO_88 * grad_0_gt0[pp];
const double DENDRO_139 = DENDRO_94 * grad_1_gt0[pp];
const double DENDRO_140 = DENDRO_20 * DENDRO_83;
const double DENDRO_141 = 2.0 * DENDRO_140;
const double DENDRO_142 = DENDRO_93 * grad_0_gt0[pp];
const double DENDRO_143 = DENDRO_94 * grad_2_gt0[pp];
const double DENDRO_144 = DENDRO_99 * grad_1_gt0[pp];
const double DENDRO_145 = DENDRO_54 * grad_0_gt3[pp];
const double DENDRO_146 = DENDRO_144 + DENDRO_145;
const double DENDRO_147 = DENDRO_108 * grad_2_gt0[pp];
const double DENDRO_148 = DENDRO_60 * grad_0_gt5[pp];
const double DENDRO_149 = DENDRO_100 * DENDRO_116;
const double DENDRO_150 = -grad_0_gt4[pp] + grad_1_gt2[pp] + grad_2_gt1[pp];
const double DENDRO_151 = 0.5 * DENDRO_150;
const double DENDRO_152 = DENDRO_149 + DENDRO_151 * DENDRO_99;
const double DENDRO_153 = DENDRO_109 * DENDRO_127;
const double DENDRO_154 = 0.5 * DENDRO_108;
const double DENDRO_155 = 0.25 * DENDRO_138;
const double DENDRO_156 = 4 * DENDRO_140;
const double DENDRO_157 = 0.25 * DENDRO_142;
const double DENDRO_158 = DENDRO_116 * grad_1_gt0[pp];
const double DENDRO_159 = DENDRO_54 * DENDRO_91;
const double DENDRO_160 = DENDRO_32 * DENDRO_83;
const double DENDRO_161 = 2.0 * DENDRO_160;
const double DENDRO_162 = DENDRO_161 * (DENDRO_158 + DENDRO_159);
const double DENDRO_163 = DENDRO_127 * grad_2_gt0[pp];
const double DENDRO_164 = DENDRO_60 * DENDRO_86;
const double DENDRO_165 = 0.5 * grad_0_gt3[pp];
const double DENDRO_166 = DENDRO_102 + DENDRO_165;
const double DENDRO_167 = DENDRO_116 * DENDRO_166;
const double DENDRO_168 = DENDRO_91 * DENDRO_99;
const double DENDRO_169 = 0.25 * DENDRO_168;
const double DENDRO_170 = -DENDRO_169;
const double DENDRO_171 = 0.5 * grad_0_gt5[pp];
const double DENDRO_172 = DENDRO_111 + DENDRO_171;
const double DENDRO_173 = -DENDRO_127;
const double DENDRO_174 = DENDRO_172 * DENDRO_173;
const double DENDRO_175 = -DENDRO_108;
const double DENDRO_176 = DENDRO_175 * DENDRO_86;
const double DENDRO_177 = -0.25 * DENDRO_176;
const double DENDRO_178 = 0.5 * DENDRO_48;
const double DENDRO_179 = DENDRO_116 * DENDRO_178;
const double DENDRO_180 = 4 * DENDRO_160;
const double DENDRO_181 = DENDRO_180 * (DENDRO_150 * DENDRO_54 + DENDRO_179);
const double DENDRO_182 = 0.5 * DENDRO_51;
const double DENDRO_183 = DENDRO_178 * DENDRO_99;
const double DENDRO_184 = -2 * DENDRO_166 * DENDRO_54 + DENDRO_183;
const double DENDRO_185 = -grad2_0_0_chi[pp];
const double DENDRO_186 = -DENDRO_94;
const double DENDRO_187 = DENDRO_30 * grad_0_chi[pp];
const double DENDRO_188 = DENDRO_30 * grad_1_chi[pp];
const double DENDRO_189 = -DENDRO_60;
const double DENDRO_190 = DENDRO_30 * grad_2_chi[pp];
const double DENDRO_191 = 2 * DENDRO_44;
const double DENDRO_192 = DENDRO_116 * DENDRO_32;
const double DENDRO_193 = DENDRO_25 * grad_2_gt3[pp];
const double DENDRO_194 = DENDRO_22 * grad_1_gt5[pp];
const double DENDRO_195 = DENDRO_150 * DENDRO_20;
const double DENDRO_196 = DENDRO_194 + DENDRO_195;
const double DENDRO_197 = -DENDRO_193 + DENDRO_196;
const double DENDRO_198 = 0.5 * grad_2_gt5[pp];
const double DENDRO_199 = DENDRO_198 * DENDRO_22;
const double DENDRO_200 = 0.5 * grad_1_gt5[pp];
const double DENDRO_201 = 1.0 * grad_2_gt4[pp];
const double DENDRO_202 = -DENDRO_201;
const double DENDRO_203 = DENDRO_200 + DENDRO_202;
const double DENDRO_204 =
    -DENDRO_172 * DENDRO_20 + DENDRO_199 + DENDRO_203 * DENDRO_25;
const double DENDRO_205 = DENDRO_204 * DENDRO_38;
const double DENDRO_206 = 0.5 * grad_1_gt3[pp];
const double DENDRO_207 = DENDRO_206 * DENDRO_25;
const double DENDRO_208 = DENDRO_166 * DENDRO_20;
const double DENDRO_209 = 1.0 * grad_1_gt4[pp];
const double DENDRO_210 = 0.5 * grad_2_gt3[pp];
const double DENDRO_211 = DENDRO_209 - DENDRO_210;
const double DENDRO_212 = DENDRO_207 + DENDRO_208 - DENDRO_211 * DENDRO_22;
const double DENDRO_213 = -DENDRO_212;
const double DENDRO_214 = DENDRO_213 * DENDRO_25;
const double DENDRO_215 = DENDRO_128 + DENDRO_205 + DENDRO_214;
const double DENDRO_216 =
    DENDRO_83 * (DENDRO_192 - 1.0 * DENDRO_197 * DENDRO_22 -
                 1.0 * DENDRO_20 * DENDRO_99 + DENDRO_215);
const double DENDRO_217 = 2.0 * grad_1_gt0[pp];
const double DENDRO_218 = -DENDRO_93;
const double DENDRO_219 = DENDRO_218 * DENDRO_32;
const double DENDRO_220 = DENDRO_150 * DENDRO_34 - DENDRO_20 * grad_2_gt3[pp] +
                          DENDRO_32 * grad_1_gt5[pp];
const double DENDRO_221 = -DENDRO_220;
const double DENDRO_222 = -DENDRO_88;
const double DENDRO_223 =
    -DENDRO_172 * DENDRO_34 + DENDRO_198 * DENDRO_32 + DENDRO_20 * DENDRO_203;
const double DENDRO_224 = -DENDRO_223;
const double DENDRO_225 = DENDRO_224 * DENDRO_38;
const double DENDRO_226 = DENDRO_20 * DENDRO_206;
const double DENDRO_227 =
    DENDRO_166 * DENDRO_34 - DENDRO_211 * DENDRO_32 + DENDRO_226;
const double DENDRO_228 = DENDRO_227 * DENDRO_25;
const double DENDRO_229 = DENDRO_186 * DENDRO_34;
const double DENDRO_230 = DENDRO_225 + DENDRO_228 + DENDRO_229;
const double DENDRO_231 =
    DENDRO_83 * (-1.0 * DENDRO_20 * DENDRO_222 + DENDRO_219 -
                 1.0 * DENDRO_22 * DENDRO_221 + DENDRO_230);
const double DENDRO_232 = 2.0 * grad_0_gt0[pp];
const double DENDRO_233 = DENDRO_175 * DENDRO_32;
const double DENDRO_234 = DENDRO_38 * grad_1_gt5[pp];
const double DENDRO_235 = DENDRO_150 * DENDRO_32;
const double DENDRO_236 = -DENDRO_22 * grad_2_gt3[pp] + DENDRO_234 + DENDRO_235;
const double DENDRO_237 = -DENDRO_236;
const double DENDRO_238 = DENDRO_198 * DENDRO_38;
const double DENDRO_239 = DENDRO_203 * DENDRO_22;
const double DENDRO_240 = -DENDRO_172 * DENDRO_32 + DENDRO_238 + DENDRO_239;
const double DENDRO_241 = -DENDRO_240;
const double DENDRO_242 = DENDRO_241 * DENDRO_38;
const double DENDRO_243 = DENDRO_206 * DENDRO_22;
const double DENDRO_244 =
    DENDRO_166 * DENDRO_32 - DENDRO_211 * DENDRO_38 + DENDRO_243;
const double DENDRO_245 = DENDRO_244 * DENDRO_25;
const double DENDRO_246 = DENDRO_189 * DENDRO_34;
const double DENDRO_247 = DENDRO_242 + DENDRO_245 + DENDRO_246;
const double DENDRO_248 =
    DENDRO_83 * (-1.0 * DENDRO_173 * DENDRO_20 - 1.0 * DENDRO_22 * DENDRO_237 +
                 DENDRO_233 + DENDRO_247);
const double DENDRO_249 = 2.0 * grad_2_gt0[pp];
const double DENDRO_250 = 3 * DENDRO_44;
const double DENDRO_251 = grad_0_chi[pp] * grad_1_chi[pp];
const double DENDRO_252 = 2 * DENDRO_20;
const double DENDRO_253 = DENDRO_250 * grad_2_chi[pp];
const double DENDRO_254 = 2 * DENDRO_32;
const double DENDRO_255 = 2 * DENDRO_22;
const double DENDRO_256 = pow(grad_1_chi[pp], 2);
const double DENDRO_257 = pow(grad_2_chi[pp], 2);
const double DENDRO_258 = -1.0 * DENDRO_192 + DENDRO_197 * DENDRO_22 +
                          DENDRO_20 * DENDRO_99 - DENDRO_215;
const double DENDRO_259 = DENDRO_20 * DENDRO_222 - 1.0 * DENDRO_219 +
                          DENDRO_22 * DENDRO_221 - DENDRO_230;
const double DENDRO_260 = DENDRO_173 * DENDRO_20 + DENDRO_22 * DENDRO_237 -
                          1.0 * DENDRO_233 - DENDRO_247;
const double DENDRO_261 =
    2 * DENDRO_187 * DENDRO_259 + 2 * DENDRO_188 * DENDRO_258 +
    2 * DENDRO_190 * DENDRO_260 +
    DENDRO_25 * (-DENDRO_250 * DENDRO_256 + 2 * grad2_1_1_chi[pp]) -
    DENDRO_252 * (-DENDRO_250 * DENDRO_251 + 2 * grad2_0_1_chi[pp]) +
    DENDRO_254 * (-DENDRO_253 * grad_0_chi[pp] + 2 * grad2_0_2_chi[pp]) -
    DENDRO_255 * (-DENDRO_253 * grad_1_chi[pp] + 2 * grad2_1_2_chi[pp]) +
    DENDRO_34 * (-DENDRO_250 * DENDRO_76 + 2 * grad2_0_0_chi[pp]) +
    DENDRO_38 * (-DENDRO_250 * DENDRO_257 + 2 * grad2_2_2_chi[pp]);
const double DENDRO_262 = DENDRO_44 * DENDRO_69;
const double DENDRO_263 = 3 * alpha[pp];
const double DENDRO_264 = DENDRO_44 * gt5[pp];
const double DENDRO_265 = DENDRO_204 + DENDRO_264 * DENDRO_43;
const double DENDRO_266 = 4 * grad_1_alpha[pp];
const double DENDRO_267 = DENDRO_266 * DENDRO_30;
const double DENDRO_268 = DENDRO_223 - 0.5 * DENDRO_44 * DENDRO_71 * gt5[pp];
const double DENDRO_269 = 4 * grad_0_alpha[pp];
const double DENDRO_270 = DENDRO_269 * DENDRO_30;
const double DENDRO_271 = 1.0 * grad_2_chi[pp];
const double DENDRO_272 = DENDRO_30 * gt5[pp];
const double DENDRO_273 = 0.5 * DENDRO_59;
const double DENDRO_274 = -DENDRO_172 * DENDRO_30 * DENDRO_32 +
                          DENDRO_238 * DENDRO_30 + DENDRO_239 * DENDRO_30 +
                          DENDRO_44 * (DENDRO_271 - DENDRO_272 * DENDRO_273);
const double DENDRO_275 = 4 * grad_2_alpha[pp];
const double DENDRO_276 = 4 * gt2[pp];
const double DENDRO_277 = 4 * gt4[pp];
const double DENDRO_278 = DENDRO_257 * DENDRO_75;
const double DENDRO_279 = DENDRO_79 * grad2_0_1_gt5[pp];
const double DENDRO_280 = DENDRO_81 * grad2_1_2_gt5[pp];
const double DENDRO_281 = DENDRO_30 * DENDRO_32;
const double DENDRO_282 = 4 * DENDRO_281;
const double DENDRO_283 = DENDRO_30 * DENDRO_34;
const double DENDRO_284 = 2.0 * DENDRO_283;
const double DENDRO_285 = DENDRO_25 * DENDRO_30;
const double DENDRO_286 = 2.0 * DENDRO_285;
const double DENDRO_287 = DENDRO_30 * DENDRO_38;
const double DENDRO_288 = 2.0 * DENDRO_287;
const double DENDRO_289 = DENDRO_175 * grad_0_gt5[pp];
const double DENDRO_290 = DENDRO_34 * DENDRO_83;
const double DENDRO_291 = 3.0 * DENDRO_290;
const double DENDRO_292 = DENDRO_237 * grad_1_gt5[pp];
const double DENDRO_293 = 3.0 * DENDRO_103;
const double DENDRO_294 = 6.0 * DENDRO_83;
const double DENDRO_295 = 0.25 * grad_2_gt3[pp];
const double DENDRO_296 = DENDRO_104 * DENDRO_197 * (DENDRO_209 - DENDRO_295);
const double DENDRO_297 = -DENDRO_136 + DENDRO_49;
const double DENDRO_298 = 4 * DENDRO_290;
const double DENDRO_299 = 0.75 * grad_0_gt4[pp];
const double DENDRO_300 = -DENDRO_123;
const double DENDRO_301 =
    DENDRO_116 * DENDRO_298 * (DENDRO_119 + DENDRO_299 + DENDRO_300);
const double DENDRO_302 = DENDRO_117 + DENDRO_122 + DENDRO_300;
const double DENDRO_303 = DENDRO_110 + DENDRO_171;
const double DENDRO_304 = DENDRO_129 * DENDRO_205 * (DENDRO_200 + DENDRO_201);
const double DENDRO_305 = DENDRO_109 * DENDRO_237;
const double DENDRO_306 = 0.25 * grad_1_gt5[pp];
const double DENDRO_307 = DENDRO_175 * DENDRO_306;
const double DENDRO_308 = DENDRO_197 * grad_1_gt5[pp];
const double DENDRO_309 = DENDRO_204 * grad_2_gt3[pp];
const double DENDRO_310 = DENDRO_308 + DENDRO_309;
const double DENDRO_311 = 2.0 * DENDRO_134;
const double DENDRO_312 = DENDRO_224 * grad_2_gt0[pp];
const double DENDRO_313 = DENDRO_175 * grad_2_gt5[pp];
const double DENDRO_314 = DENDRO_241 * grad_0_gt5[pp];
const double DENDRO_315 = DENDRO_237 * grad_2_gt5[pp];
const double DENDRO_316 = DENDRO_241 * grad_1_gt5[pp];
const double DENDRO_317 = DENDRO_116 * DENDRO_295;
const double DENDRO_318 = 0.5 * DENDRO_86;
const double DENDRO_319 = DENDRO_197 * DENDRO_318 + DENDRO_317;
const double DENDRO_320 = 0.25 * DENDRO_313;
const double DENDRO_321 = 0.25 * DENDRO_315;
const double DENDRO_322 = DENDRO_136 * DENDRO_221;
const double DENDRO_323 = DENDRO_116 * grad_1_gt5[pp];
const double DENDRO_324 = DENDRO_204 * DENDRO_91;
const double DENDRO_325 = DENDRO_161 * (DENDRO_323 + DENDRO_324);
const double DENDRO_326 = DENDRO_116 * DENDRO_211;
const double DENDRO_327 = DENDRO_197 * DENDRO_91;
const double DENDRO_328 = 0.25 * DENDRO_327;
const double DENDRO_329 = DENDRO_326 + DENDRO_328;
const double DENDRO_330 = DENDRO_150 * DENDRO_224;
const double DENDRO_331 = DENDRO_221 * DENDRO_51;
const double DENDRO_332 = DENDRO_150 * DENDRO_218;
const double DENDRO_333 = 0.25 * DENDRO_332;
const double DENDRO_334 = 0.5 * DENDRO_203;
const double DENDRO_335 = DENDRO_116 * DENDRO_334;
const double DENDRO_336 = 0.5 * DENDRO_172;
const double DENDRO_337 = DENDRO_221 * DENDRO_336;
const double DENDRO_338 = -DENDRO_197 * DENDRO_334;
const double DENDRO_339 = 2 * DENDRO_204 * DENDRO_211 + DENDRO_338;
const double DENDRO_340 = -DENDRO_218 * DENDRO_336;
const double DENDRO_341 = 2 * DENDRO_51;
const double DENDRO_342 = -grad2_2_2_chi[pp];
const double DENDRO_343 = -DENDRO_32;
const double DENDRO_344 = -DENDRO_172;
const double DENDRO_345 = -DENDRO_34;
const double DENDRO_346 = -DENDRO_203;
const double DENDRO_347 =
    DENDRO_198 * DENDRO_343 + DENDRO_20 * DENDRO_346 + DENDRO_344 * DENDRO_345;
const double DENDRO_348 = -DENDRO_25;
const double DENDRO_349 =
    DENDRO_199 + DENDRO_20 * DENDRO_344 + DENDRO_346 * DENDRO_348;
const double DENDRO_350 = -DENDRO_38;
const double DENDRO_351 = DENDRO_198 * DENDRO_350;
const double DENDRO_352 = DENDRO_343 * DENDRO_344;
const double DENDRO_353 = DENDRO_22 * DENDRO_346;
const double DENDRO_354 = DENDRO_258 * DENDRO_83;
const double DENDRO_355 = 2.0 * DENDRO_354;
const double DENDRO_356 = DENDRO_259 * DENDRO_83;
const double DENDRO_357 = 2.0 * DENDRO_356;
const double DENDRO_358 = DENDRO_260 * DENDRO_83;
const double DENDRO_359 = 2.0 * DENDRO_358;
const double DENDRO_360 = -DENDRO_261;
const double DENDRO_361 = DENDRO_360 * DENDRO_44;
const double DENDRO_362 = DENDRO_44 * gt3[pp];
const double DENDRO_363 = DENDRO_275 * DENDRO_30;
const double DENDRO_364 = 1.0 * grad_1_chi[pp];
const double DENDRO_365 = DENDRO_30 * gt3[pp];
const double DENDRO_366 = DENDRO_207 * DENDRO_30 + DENDRO_208 * DENDRO_30 -
                          DENDRO_211 * DENDRO_22 * DENDRO_30 +
                          DENDRO_44 * (DENDRO_364 - DENDRO_365 * DENDRO_43);
const double DENDRO_367 = -grad2_1_1_chi[pp];
const double DENDRO_368 = -DENDRO_166;
const double DENDRO_369 =
    DENDRO_211 * DENDRO_343 + DENDRO_226 + DENDRO_345 * DENDRO_368;
const double DENDRO_370 = DENDRO_206 * DENDRO_348;
const double DENDRO_371 = DENDRO_20 * DENDRO_368;
const double DENDRO_372 = DENDRO_211 * DENDRO_22;
const double DENDRO_373 =
    DENDRO_211 * DENDRO_350 + DENDRO_243 + DENDRO_343 * DENDRO_368;
const double DENDRO_374 = DENDRO_221 * DENDRO_48;
const double DENDRO_375 = DENDRO_150 * DENDRO_222;
const double DENDRO_376 = 0.25 * DENDRO_375;
const double DENDRO_377 = DENDRO_132 * DENDRO_221;
const double DENDRO_378 = 0.5 * DENDRO_91;
const double DENDRO_379 = DENDRO_173 * DENDRO_306;
const double DENDRO_380 = DENDRO_173 * DENDRO_203;
const double DENDRO_381 = DENDRO_237 * DENDRO_86;
const double DENDRO_382 = -0.25 * DENDRO_381;
const double DENDRO_383 = 0.5 * DENDRO_166;
const double DENDRO_384 = DENDRO_221 * DENDRO_383;
const double DENDRO_385 = 0.5 * DENDRO_211;
const double DENDRO_386 = DENDRO_237 * DENDRO_385;
const double DENDRO_387 = 2 * DENDRO_203 * DENDRO_244;
const double DENDRO_388 = DENDRO_197 * grad_1_gt3[pp];
const double DENDRO_389 = 0.25 * DENDRO_388;
const double DENDRO_390 = DENDRO_213 * grad_2_gt3[pp];
const double DENDRO_391 = DENDRO_173 * DENDRO_385;
const double DENDRO_392 = -DENDRO_222 * DENDRO_383;
const double DENDRO_393 = 2 * DENDRO_227 * DENDRO_48;
const double DENDRO_394 = DENDRO_99 * grad_1_gt3[pp];
const double DENDRO_395 = 0.25 * DENDRO_394;
const double DENDRO_396 = DENDRO_213 * grad_0_gt3[pp];
const double DENDRO_397 = DENDRO_227 * grad_1_gt0[pp];
const double DENDRO_398 = DENDRO_150 * DENDRO_227;
const double DENDRO_399 = DENDRO_244 * grad_1_gt5[pp];
const double DENDRO_400 = DENDRO_244 * DENDRO_86;
const double DENDRO_401 = DENDRO_202 + DENDRO_306;
const double DENDRO_402 = -DENDRO_119;
const double DENDRO_403 = DENDRO_114 * (DENDRO_117 + DENDRO_120 + DENDRO_402);
const double DENDRO_404 = DENDRO_298 * (-DENDRO_132 + DENDRO_46);
const double DENDRO_405 = DENDRO_298 * (DENDRO_123 + DENDRO_299 + DENDRO_402);
const double DENDRO_406 = 4 * gt1[pp];
const double DENDRO_407 = DENDRO_100 * DENDRO_197;
const double DENDRO_408 = DENDRO_295 * DENDRO_99;
const double DENDRO_409 = DENDRO_99 * grad_0_gt3[pp];
const double DENDRO_410 = DENDRO_197 * grad_2_gt3[pp];
const double DENDRO_411 = 3.0 * DENDRO_113;
const double DENDRO_412 =
    -DENDRO_129 * DENDRO_228 * (DENDRO_101 + DENDRO_165) -
    DENDRO_129 * DENDRO_245 * (DENDRO_209 + DENDRO_210) -
    DENDRO_180 * (DENDRO_165 * DENDRO_197 + DENDRO_408) -
    DENDRO_180 * (DENDRO_210 * DENDRO_99 + DENDRO_407) -
    DENDRO_256 * DENDRO_75 + DENDRO_277 * grad_1_Gt2[pp] +
    DENDRO_282 * grad2_0_2_gt3[pp] + DENDRO_284 * grad2_0_0_gt3[pp] +
    DENDRO_286 * grad2_1_1_gt3[pp] + DENDRO_288 * grad2_2_2_gt3[pp] -
    DENDRO_291 * DENDRO_409 + DENDRO_406 * grad_1_Gt0[pp] -
    DENDRO_410 * DENDRO_411 - DENDRO_79 * grad2_0_1_gt3[pp] -
    DENDRO_81 * grad2_1_2_gt3[pp] + 4 * grad_1_Gt1[pp] * gt3[pp];
const double DENDRO_413 = DENDRO_222 * grad_1_gt0[pp];
const double DENDRO_414 = DENDRO_218 * grad_2_gt0[pp];
const double DENDRO_415 = DENDRO_132 * DENDRO_218;
const double DENDRO_416 = DENDRO_136 * DENDRO_222;
const double DENDRO_417 = DENDRO_222 * grad_0_gt0[pp];
const double DENDRO_418 = DENDRO_186 * grad_1_gt0[pp];
const double DENDRO_419 = DENDRO_218 * grad_0_gt0[pp];
const double DENDRO_420 = DENDRO_186 * grad_2_gt0[pp];
const double DENDRO_421 = DENDRO_189 * grad_0_gt5[pp];
const double DENDRO_422 = 0.25 * DENDRO_417;
const double DENDRO_423 = 0.25 * DENDRO_419;
const double DENDRO_424 = DENDRO_109 * DENDRO_173;
const double DENDRO_425 = DENDRO_189 * DENDRO_86;
const double DENDRO_426 = DENDRO_173 * DENDRO_182;
const double DENDRO_427 = DENDRO_175 * DENDRO_182;
const double DENDRO_428 = DENDRO_345 * DENDRO_52;
const double DENDRO_429 = DENDRO_343 * DENDRO_51;
const double DENDRO_430 = DENDRO_348 * DENDRO_48 + DENDRO_53;
const double DENDRO_431 =
    DENDRO_22 * DENDRO_48 + DENDRO_343 * DENDRO_52 + DENDRO_350 * DENDRO_51;
const double DENDRO_432 = DENDRO_197 * grad_1_gt0[pp];
const double DENDRO_433 = DENDRO_116 * grad_0_gt3[pp];
const double DENDRO_434 = DENDRO_86 * DENDRO_99;
const double DENDRO_435 = DENDRO_433 + DENDRO_434;
const double DENDRO_436 = DENDRO_237 * grad_2_gt0[pp];
const double DENDRO_437 = DENDRO_173 * grad_0_gt5[pp];
const double DENDRO_438 = DENDRO_176 + DENDRO_437;
const double DENDRO_439 = DENDRO_109 * DENDRO_218;
const double DENDRO_440 = -DENDRO_172 * DENDRO_241;
const double DENDRO_441 = DENDRO_151 * DENDRO_204 + 0.25 * DENDRO_323;
const double DENDRO_442 = DENDRO_222 * DENDRO_86;
const double DENDRO_443 = 0.25 * DENDRO_442;
const double DENDRO_444 = DENDRO_385 * DENDRO_99;
const double DENDRO_445 = -DENDRO_197 * DENDRO_383;
const double DENDRO_446 = DENDRO_444 + DENDRO_445;
const double DENDRO_447 = DENDRO_186 * DENDRO_51;
const double DENDRO_448 = 0.25 * DENDRO_158 + DENDRO_318 * DENDRO_54;
const double DENDRO_449 = DENDRO_136 * DENDRO_175;
const double DENDRO_450 = 0.25 * DENDRO_116;
const double DENDRO_451 = DENDRO_150 * DENDRO_450 + DENDRO_200 * DENDRO_54;
const double DENDRO_452 = DENDRO_189 * DENDRO_198;
const double DENDRO_453 = -DENDRO_175 * DENDRO_336;
const double DENDRO_454 = DENDRO_136 * DENDRO_218;
const double DENDRO_455 = DENDRO_224 * DENDRO_52;
const double DENDRO_456 = DENDRO_182 * DENDRO_218 + DENDRO_455;
const double DENDRO_457 = DENDRO_450 * DENDRO_91;
const double DENDRO_458 = DENDRO_204 * DENDRO_47 + DENDRO_450 * DENDRO_86;
const double DENDRO_459 = DENDRO_166 * DENDRO_204;
const double DENDRO_460 = DENDRO_224 * DENDRO_47;
const double DENDRO_461 = DENDRO_109 * DENDRO_222 + DENDRO_460;
const double DENDRO_462 = DENDRO_151 * DENDRO_241;
const double DENDRO_463 = DENDRO_305 + DENDRO_307;
const double DENDRO_464 = DENDRO_150 * DENDRO_197;
const double DENDRO_465 = 0.25 * DENDRO_464;
const double DENDRO_466 = DENDRO_165 * DENDRO_204;
const double DENDRO_467 = DENDRO_306 * DENDRO_99 + DENDRO_466;
const double DENDRO_468 = 0.25 * DENDRO_86;
const double DENDRO_469 = DENDRO_218 * DENDRO_468 + DENDRO_460;
const double DENDRO_470 = DENDRO_237 * DENDRO_336;
const double DENDRO_471 = -DENDRO_470;
const double DENDRO_472 = 0.25 * DENDRO_173;
const double DENDRO_473 = DENDRO_472 * grad_2_gt5[pp];
const double DENDRO_474 = DENDRO_241 * DENDRO_318 + DENDRO_473;
const double DENDRO_475 = DENDRO_178 * DENDRO_197;
const double DENDRO_476 = -0.5 * DENDRO_167 + DENDRO_211 * DENDRO_54;
const double DENDRO_477 = 0.25 * DENDRO_221;
const double DENDRO_478 = DENDRO_477 * grad_0_gt0[pp];
const double DENDRO_479 = DENDRO_182 * DENDRO_222;
const double DENDRO_480 = DENDRO_478 + DENDRO_479;
const double DENDRO_481 = DENDRO_186 * DENDRO_318 + DENDRO_415;
const double DENDRO_482 = DENDRO_189 * DENDRO_200;
const double DENDRO_483 = 0.25 * DENDRO_150;
const double DENDRO_484 = DENDRO_175 * DENDRO_483;
const double DENDRO_485 = DENDRO_150 * DENDRO_237;
const double DENDRO_486 = DENDRO_173 * grad_1_gt5[pp];
const double DENDRO_487 = DENDRO_381 + DENDRO_486;
const double DENDRO_488 = 1.0 * DENDRO_103;
const double DENDRO_489 = -grad2_0_2_chi[pp];
const double DENDRO_490 = DENDRO_343 * grad_0_gt5[pp];
const double DENDRO_491 = DENDRO_345 * grad_2_gt0[pp];
const double DENDRO_492 = 0.5 * DENDRO_187;
const double DENDRO_493 = DENDRO_115 + DENDRO_348 * DENDRO_91;
const double DENDRO_494 = 0.5 * DENDRO_188;
const double DENDRO_495 = DENDRO_350 * grad_0_gt5[pp];
const double DENDRO_496 = DENDRO_343 * grad_2_gt0[pp];
const double DENDRO_497 = DENDRO_22 * DENDRO_91;
const double DENDRO_498 = 0.5 * DENDRO_190;
const double DENDRO_499 = 2.0 * gt2[pp];
const double DENDRO_500 = 2.0 * gt4[pp];
const double DENDRO_501 = 2.0 * gt5[pp];
const double DENDRO_502 = 2.0 * grad_2_Gt0[pp];
const double DENDRO_503 = 2.0 * grad_2_Gt1[pp];
const double DENDRO_504 = DENDRO_75 * grad_2_chi[pp];
const double DENDRO_505 = DENDRO_504 * grad_0_chi[pp];
const double DENDRO_506 = DENDRO_79 * grad2_0_1_gt2[pp];
const double DENDRO_507 = DENDRO_81 * grad2_1_2_gt2[pp];
const double DENDRO_508 = DENDRO_30 * gt2[pp];
const double DENDRO_509 =
    -DENDRO_191 *
        (DENDRO_489 + DENDRO_492 * (DENDRO_490 + DENDRO_491 + DENDRO_92) +
         DENDRO_493 * DENDRO_494 +
         DENDRO_498 * (DENDRO_495 + DENDRO_496 + DENDRO_497)) +
    DENDRO_282 * grad2_0_2_gt2[pp] + DENDRO_284 * grad2_0_0_gt2[pp] +
    DENDRO_286 * grad2_1_1_gt2[pp] + DENDRO_288 * grad2_2_2_gt2[pp] +
    DENDRO_355 * grad_1_gt2[pp] + DENDRO_357 * grad_0_gt2[pp] +
    DENDRO_359 * grad_2_gt2[pp] + DENDRO_361 * DENDRO_508 +
    DENDRO_499 * grad_0_Gt0[pp] + DENDRO_499 * grad_2_Gt2[pp] +
    DENDRO_500 * grad_0_Gt1[pp] + DENDRO_501 * grad_0_Gt2[pp] +
    DENDRO_502 * gt0[pp] + DENDRO_503 * gt1[pp] - DENDRO_505 - DENDRO_506 -
    DENDRO_507;
const double DENDRO_510 =
    -DENDRO_20 * DENDRO_30 * DENDRO_91 + DENDRO_30 * DENDRO_89 +
    DENDRO_30 * DENDRO_90 +
    DENDRO_44 * (-DENDRO_508 * DENDRO_71 + grad_2_chi[pp]);
const double DENDRO_511 = 2.0 * grad_0_alpha[pp];
const double DENDRO_512 =
    DENDRO_106 * DENDRO_30 + DENDRO_107 * DENDRO_30 -
    DENDRO_22 * DENDRO_30 * DENDRO_91 +
    DENDRO_44 * (-DENDRO_508 * DENDRO_59 + grad_0_chi[pp]);
const double DENDRO_513 = 2.0 * grad_2_alpha[pp];
const double DENDRO_514 = 2.0 * grad_1_alpha[pp];
const double DENDRO_515 = DENDRO_44 * gt2[pp];
const double DENDRO_516 = DENDRO_30 * (DENDRO_116 + DENDRO_42 * DENDRO_515);
const double DENDRO_517 = -DENDRO_510 * DENDRO_511 - DENDRO_512 * DENDRO_513 +
                          DENDRO_514 * DENDRO_516 - 4 * grad2_0_2_alpha[pp];
const double DENDRO_518 = DENDRO_221 * grad_2_gt0[pp];
const double DENDRO_519 = DENDRO_222 * grad_0_gt5[pp] + DENDRO_518;
const double DENDRO_520 = DENDRO_116 * grad_2_gt3[pp];
const double DENDRO_521 = DENDRO_520 + DENDRO_99 * grad_1_gt5[pp];
const double DENDRO_522 = DENDRO_464 + DENDRO_521;
const double DENDRO_523 = DENDRO_114 * (-DENDRO_335 + DENDRO_441);
const double DENDRO_524 = 0.25 * DENDRO_485;
const double DENDRO_525 = DENDRO_104 * (DENDRO_408 + DENDRO_446);
const double DENDRO_526 = DENDRO_298 * (0.5 * DENDRO_159 + DENDRO_448);
const double DENDRO_527 = DENDRO_180 * (-DENDRO_203 * DENDRO_54 + DENDRO_458);
const double DENDRO_528 = DENDRO_109 * DENDRO_175 + DENDRO_452;
const double DENDRO_529 = DENDRO_180 * (DENDRO_451 + DENDRO_457);
const double DENDRO_530 = DENDRO_334 * DENDRO_99;
const double DENDRO_531 =
    -0.5 * DENDRO_116 * DENDRO_211 + DENDRO_459 + DENDRO_530;
const double DENDRO_532 = DENDRO_307 + DENDRO_473;
const double DENDRO_533 = -DENDRO_222 * DENDRO_336;
const double DENDRO_534 = DENDRO_169 + DENDRO_476;
const double DENDRO_535 = DENDRO_151 * DENDRO_186;
const double DENDRO_536 = DENDRO_136 * DENDRO_237;
const double DENDRO_537 = DENDRO_424 + DENDRO_482;
const double DENDRO_538 = 0.25 * DENDRO_434;
const double DENDRO_539 = DENDRO_210 * DENDRO_54;
const double DENDRO_540 = DENDRO_132 * DENDRO_197 + DENDRO_539;
const double DENDRO_541 = DENDRO_538 + DENDRO_540;
const double DENDRO_542 = DENDRO_221 * grad_1_gt0[pp];
const double DENDRO_543 = DENDRO_442 + DENDRO_542;
const double DENDRO_544 = DENDRO_218 * grad_0_gt3[pp];
const double DENDRO_545 = DENDRO_175 * grad_2_gt3[pp];
const double DENDRO_546 = 0.25 * DENDRO_308;
const double DENDRO_547 = -DENDRO_203 * DENDRO_241;
const double DENDRO_548 = DENDRO_109 * DENDRO_221 + DENDRO_224 * DENDRO_378;
const double DENDRO_549 = DENDRO_211 * DENDRO_213;
const double DENDRO_550 = DENDRO_227 * DENDRO_318;
const double DENDRO_551 = DENDRO_100 * DENDRO_221 + DENDRO_550;
const double DENDRO_552 = DENDRO_237 * DENDRO_295;
const double DENDRO_553 = DENDRO_178 * DENDRO_218;
const double DENDRO_554 = DENDRO_479 + DENDRO_553;
const double DENDRO_555 = DENDRO_224 * DENDRO_48 + 0.5 * DENDRO_331;
const double DENDRO_556 = DENDRO_218 * DENDRO_91;
const double DENDRO_557 = 0.25 * DENDRO_556;
const double DENDRO_558 = DENDRO_241 * DENDRO_378;
const double DENDRO_559 = DENDRO_197 * DENDRO_468 + DENDRO_466;
const double DENDRO_560 = DENDRO_175 * DENDRO_334;
const double DENDRO_561 = -DENDRO_560;
const double DENDRO_562 = DENDRO_198 * DENDRO_244;
const double DENDRO_563 = -DENDRO_237 * DENDRO_334 + DENDRO_562;
const double DENDRO_564 = DENDRO_171 * DENDRO_227 + DENDRO_477 * DENDRO_91;
const double DENDRO_565 = DENDRO_204 * DENDRO_206;
const double DENDRO_566 = DENDRO_197 * DENDRO_295 + DENDRO_565;
const double DENDRO_567 = DENDRO_197 * DENDRO_385;
const double DENDRO_568 = DENDRO_150 * DENDRO_477;
const double DENDRO_569 = DENDRO_165 * DENDRO_224 + DENDRO_221 * DENDRO_468;
const double DENDRO_570 = -DENDRO_218 * DENDRO_383;
const double DENDRO_571 = DENDRO_227 * DENDRO_51 + 0.5 * DENDRO_374;
const double DENDRO_572 = DENDRO_450 * grad_1_gt3[pp];
const double DENDRO_573 = DENDRO_407 + DENDRO_572;
const double DENDRO_574 = DENDRO_213 * DENDRO_318;
const double DENDRO_575 = 0.25 * DENDRO_91;
const double DENDRO_576 = DENDRO_171 * DENDRO_244;
const double DENDRO_577 = DENDRO_237 * DENDRO_575 + DENDRO_576;
const double DENDRO_578 = DENDRO_175 * DENDRO_91;
const double DENDRO_579 = 1.0 * DENDRO_290;
const double DENDRO_580 = -grad2_1_2_chi[pp];
const double DENDRO_581 = DENDRO_150 * DENDRO_345 + DENDRO_20 * grad_2_gt3[pp] +
                          DENDRO_343 * grad_1_gt5[pp];
const double DENDRO_582 = DENDRO_348 * grad_2_gt3[pp];
const double DENDRO_583 = DENDRO_350 * grad_1_gt5[pp];
const double DENDRO_584 = DENDRO_22 * grad_2_gt3[pp];
const double DENDRO_585 = DENDRO_150 * DENDRO_343;
const double DENDRO_586 = DENDRO_30 * gt4[pp];
const double DENDRO_587 =
    DENDRO_282 * grad2_0_2_gt4[pp] + DENDRO_284 * grad2_0_0_gt4[pp] +
    DENDRO_286 * grad2_1_1_gt4[pp] + DENDRO_288 * grad2_2_2_gt4[pp] +
    DENDRO_499 * grad_1_Gt0[pp] + DENDRO_500 * grad_1_Gt1[pp] +
    DENDRO_500 * grad_2_Gt2[pp] + DENDRO_501 * grad_1_Gt2[pp] +
    DENDRO_502 * gt1[pp] + DENDRO_503 * gt3[pp] - DENDRO_504 * grad_1_chi[pp] -
    DENDRO_79 * grad2_0_1_gt4[pp] - DENDRO_81 * grad2_1_2_gt4[pp];
const double DENDRO_588 =
    -DENDRO_191 *
        (DENDRO_492 * DENDRO_581 + DENDRO_494 * (DENDRO_196 + DENDRO_582) +
         DENDRO_498 * (DENDRO_583 + DENDRO_584 + DENDRO_585) + DENDRO_580) +
    DENDRO_355 * grad_1_gt4[pp] + DENDRO_357 * grad_0_gt4[pp] +
    DENDRO_359 * grad_2_gt4[pp] + DENDRO_361 * DENDRO_586 + DENDRO_587;
const double DENDRO_589 = DENDRO_194 * DENDRO_30 + DENDRO_195 * DENDRO_30;
const double DENDRO_590 =
    -DENDRO_193 * DENDRO_30 -
    DENDRO_44 * (-DENDRO_42 * DENDRO_586 + grad_2_chi[pp]) + DENDRO_589;
const double DENDRO_591 =
    -DENDRO_22 * DENDRO_30 * grad_2_gt3[pp] + DENDRO_234 * DENDRO_30 +
    DENDRO_235 * DENDRO_30 +
    DENDRO_44 * (-DENDRO_586 * DENDRO_59 + grad_1_chi[pp]);
const double DENDRO_592 = DENDRO_220 - DENDRO_44 * DENDRO_71 * gt4[pp];
const double DENDRO_593 = -DENDRO_30 * DENDRO_511 * DENDRO_592 -
                          DENDRO_513 * DENDRO_591 + DENDRO_514 * DENDRO_590 -
                          4 * grad2_1_2_alpha[pp];
const double DENDRO_594 = 1.0 * DENDRO_399;
const double DENDRO_595 = 0.5 * DENDRO_398;
const double DENDRO_596 = 0.25 * DENDRO_578;
const double DENDRO_597 = DENDRO_305 + DENDRO_473;
const double DENDRO_598 = DENDRO_172 * DENDRO_227;
const double DENDRO_599 = DENDRO_565 + DENDRO_567;
const double DENDRO_600 = DENDRO_237 * DENDRO_306;
const double DENDRO_601 = DENDRO_407 + DENDRO_408;
const double DENDRO_602 = DENDRO_227 * DENDRO_50;
const double DENDRO_603 = DENDRO_100 * DENDRO_218 + DENDRO_602;
const double DENDRO_604 = DENDRO_213 * DENDRO_378 + DENDRO_572;
const double DENDRO_605 = DENDRO_175 * DENDRO_295 + DENDRO_576;
const double DENDRO_606 = 1.0 * DENDRO_160;
const double DENDRO_607 =
    -DENDRO_114 * (DENDRO_204 * DENDRO_210 + DENDRO_338 + DENDRO_546) -
    DENDRO_180 * (-DENDRO_530 + DENDRO_559) -
    DENDRO_579 * (DENDRO_168 + DENDRO_435) -
    DENDRO_606 * (DENDRO_327 + DENDRO_521);
const double DENDRO_608 = DENDRO_470 + DENDRO_560;
const double DENDRO_609 = DENDRO_100 * DENDRO_222;
const double DENDRO_610 = -DENDRO_166 * DENDRO_213;
const double DENDRO_611 = DENDRO_151 * DENDRO_244 + DENDRO_173 * DENDRO_295;
const double DENDRO_612 = DENDRO_186 * DENDRO_48;
const double DENDRO_613 = 0.25 * DENDRO_144;
const double DENDRO_614 = DENDRO_136 * DENDRO_173 + DENDRO_189 * DENDRO_378;
const double DENDRO_615 = 0.5 * DENDRO_174 + DENDRO_189 * DENDRO_203;
const double DENDRO_616 = DENDRO_478 + DENDRO_553;
const double DENDRO_617 = DENDRO_186 * DENDRO_378 + DENDRO_416;
const double DENDRO_618 = DENDRO_483 * DENDRO_99 + DENDRO_539;
const double DENDRO_619 = DENDRO_172 * DENDRO_244 + 0.5 * DENDRO_380;
const double DENDRO_620 = DENDRO_151 * DENDRO_213;
const double DENDRO_621 = DENDRO_222 * DENDRO_575 + DENDRO_602;
const double DENDRO_622 = -DENDRO_383 * DENDRO_99;
const double DENDRO_623 = DENDRO_206 * DENDRO_54;
const double DENDRO_624 = DENDRO_150 * DENDRO_472 + DENDRO_189 * DENDRO_210;
const double DENDRO_625 = DENDRO_132 * DENDRO_222;
const double DENDRO_626 = DENDRO_227 * DENDRO_52;
const double DENDRO_627 = DENDRO_178 * DENDRO_222 + DENDRO_626;
const double DENDRO_628 = DENDRO_472 * DENDRO_86;
const double DENDRO_629 = DENDRO_244 * DENDRO_50 + DENDRO_472 * DENDRO_91;
const double DENDRO_630 = 1.0 * DENDRO_113;
const double DENDRO_631 = -grad2_0_1_chi[pp];
const double DENDRO_632 = DENDRO_345 * grad_1_gt0[pp];
const double DENDRO_633 = DENDRO_343 * DENDRO_86;
const double DENDRO_634 = DENDRO_348 * grad_0_gt3[pp];
const double DENDRO_635 = DENDRO_22 * grad_0_gt3[pp];
const double DENDRO_636 =
    DENDRO_343 * grad_1_gt0[pp] + DENDRO_350 * DENDRO_86 + DENDRO_635;
const double DENDRO_637 = 2.0 * gt1[pp];
const double DENDRO_638 = DENDRO_637 * grad_0_Gt0[pp];
const double DENDRO_639 = 2.0 * grad_0_Gt1[pp] * gt3[pp];
const double DENDRO_640 = DENDRO_500 * grad_0_Gt2[pp];
const double DENDRO_641 = 2.0 * grad_1_Gt0[pp] * gt0[pp];
const double DENDRO_642 = DENDRO_637 * grad_1_Gt1[pp];
const double DENDRO_643 = DENDRO_499 * grad_1_Gt2[pp];
const double DENDRO_644 = -DENDRO_251 * DENDRO_75;
const double DENDRO_645 = DENDRO_282 * grad2_0_2_gt1[pp];
const double DENDRO_646 = DENDRO_284 * grad2_0_0_gt1[pp];
const double DENDRO_647 = DENDRO_286 * grad2_1_1_gt1[pp];
const double DENDRO_648 = DENDRO_288 * grad2_2_2_gt1[pp];
const double DENDRO_649 = -DENDRO_79 * grad2_0_1_gt1[pp];
const double DENDRO_650 = -DENDRO_81 * grad2_1_2_gt1[pp];
const double DENDRO_651 = DENDRO_30 * gt1[pp];
const double DENDRO_652 =
    -DENDRO_191 * (DENDRO_492 * (DENDRO_632 + DENDRO_633 + DENDRO_84) +
                   DENDRO_494 * (DENDRO_634 + DENDRO_98) +
                   DENDRO_498 * DENDRO_636 + DENDRO_631) +
    DENDRO_355 * grad_1_gt1[pp] + DENDRO_357 * grad_0_gt1[pp] +
    DENDRO_359 * grad_2_gt1[pp] + DENDRO_361 * DENDRO_651 + DENDRO_638 +
    DENDRO_639 + DENDRO_640 + DENDRO_641 + DENDRO_642 + DENDRO_643 +
    DENDRO_644 + DENDRO_645 + DENDRO_646 + DENDRO_647 + DENDRO_648 +
    DENDRO_649 + DENDRO_650;
const double DENDRO_653 =
    -DENDRO_20 * DENDRO_30 * grad_0_gt3[pp] + DENDRO_30 * DENDRO_85 +
    DENDRO_30 * DENDRO_87 +
    DENDRO_44 * (-DENDRO_651 * DENDRO_71 + grad_1_chi[pp]);
const double DENDRO_654 =
    -DENDRO_20 * DENDRO_30 * grad_1_gt0[pp] -
    DENDRO_22 * DENDRO_30 * DENDRO_86 + DENDRO_30 * DENDRO_95 +
    DENDRO_44 * (-DENDRO_42 * DENDRO_651 + grad_0_chi[pp]);
const double DENDRO_655 = DENDRO_44 * gt1[pp];
const double DENDRO_656 =
    DENDRO_30 * DENDRO_513 *
        (-DENDRO_125 - DENDRO_126 + DENDRO_59 * DENDRO_655 + DENDRO_635) -
    DENDRO_511 * DENDRO_653 - DENDRO_514 * DENDRO_654 - 4 * grad2_0_1_alpha[pp];
const double DENDRO_657 = -0.25 * DENDRO_175 * grad_1_gt5[pp] + DENDRO_608;
const double DENDRO_658 = 0.5 * DENDRO_394;
const double DENDRO_659 = DENDRO_227 * DENDRO_47;
const double DENDRO_660 = DENDRO_177 + DENDRO_615;
const double DENDRO_661 = -0.5 * DENDRO_175 * DENDRO_211 + DENDRO_619;
const double DENDRO_662 = DENDRO_408 + DENDRO_572;
const double DENDRO_663 = DENDRO_100 * DENDRO_99 + DENDRO_623;
const double DENDRO_664 = -DENDRO_114 * (DENDRO_116 * DENDRO_210 + DENDRO_465) +
                          DENDRO_135 * (DENDRO_445 + DENDRO_662) +
                          DENDRO_156 * (DENDRO_622 + DENDRO_663) -
                          DENDRO_180 * (DENDRO_149 + DENDRO_540) -
                          DENDRO_180 * (DENDRO_149 + DENDRO_618) -
                          DENDRO_298 * (1.0 * DENDRO_145 + DENDRO_613);
const double DENDRO_665 =
    -DENDRO_20 *
        (DENDRO_656 +
         alpha[pp] *
             (-DENDRO_104 * (DENDRO_391 + DENDRO_611) -
              DENDRO_104 * (DENDRO_610 + DENDRO_658) -
              DENDRO_104 * (DENDRO_392 + DENDRO_609 + DENDRO_659) +
              DENDRO_114 * DENDRO_657 +
              DENDRO_134 * (DENDRO_375 + DENDRO_542 + DENDRO_544) +
              DENDRO_134 * (DENDRO_485 + DENDRO_486 + DENDRO_545) -
              DENDRO_135 * DENDRO_661 + DENDRO_135 * (DENDRO_570 + DENDRO_621) +
              DENDRO_135 * (DENDRO_620 + DENDRO_662) +
              DENDRO_141 * (DENDRO_186 * grad_0_gt3[pp] + DENDRO_413) +
              DENDRO_156 * (DENDRO_624 + DENDRO_628) +
              DENDRO_156 * (-DENDRO_166 * DENDRO_186 + DENDRO_627) +
              DENDRO_156 * (DENDRO_189 * DENDRO_211 + DENDRO_629) +
              DENDRO_156 * (DENDRO_213 * DENDRO_47 + DENDRO_663) +
              DENDRO_180 * DENDRO_660 - DENDRO_180 * (DENDRO_415 + DENDRO_617) -
              DENDRO_180 * (DENDRO_535 + DENDRO_616) -
              DENDRO_180 * (DENDRO_482 + DENDRO_536 + DENDRO_596) -
              DENDRO_298 * (0.5 * DENDRO_425 + DENDRO_614) -
              DENDRO_298 * (DENDRO_186 * DENDRO_47 + DENDRO_422 + DENDRO_612) -
              DENDRO_630 * (DENDRO_332 + DENDRO_518 + DENDRO_556) + DENDRO_652 +
              DENDRO_664)) -
    DENDRO_20 *
        (DENDRO_656 +
         alpha[pp] *
             (-DENDRO_104 * (1.0 * DENDRO_397 + DENDRO_609) -
              DENDRO_104 * (0.5 * DENDRO_400 + DENDRO_611) -
              DENDRO_104 * (DENDRO_165 * DENDRO_213 + DENDRO_395 + DENDRO_610) -
              DENDRO_114 * (DENDRO_221 * DENDRO_50 + DENDRO_557) -
              DENDRO_114 * (0.25 * DENDRO_237 * grad_0_gt5[pp] - DENDRO_608) +
              DENDRO_135 * (DENDRO_377 + DENDRO_603) +
              DENDRO_135 * (DENDRO_377 + DENDRO_621) +
              DENDRO_135 * (-DENDRO_382 - DENDRO_619) +
              DENDRO_135 * (DENDRO_445 + DENDRO_604) +
              DENDRO_135 * (DENDRO_524 + DENDRO_605) +
              DENDRO_135 * (DENDRO_601 + DENDRO_620) +
              DENDRO_141 * (DENDRO_213 * grad_1_gt0[pp] + DENDRO_409) +
              DENDRO_156 * (DENDRO_625 + DENDRO_627) +
              DENDRO_156 * (DENDRO_628 + DENDRO_629) +
              DENDRO_156 * (DENDRO_244 * DENDRO_51 + DENDRO_624) +
              DENDRO_156 * (DENDRO_165 * DENDRO_186 + DENDRO_625 + DENDRO_626) +
              DENDRO_156 * (DENDRO_213 * DENDRO_48 + DENDRO_622 + DENDRO_623) -
              DENDRO_180 * (DENDRO_416 + DENDRO_616) -
              DENDRO_180 * (DENDRO_475 + DENDRO_618) -
              DENDRO_180 * (DENDRO_478 + DENDRO_617) -
              DENDRO_180 * (0.5 * DENDRO_237 * DENDRO_51 - DENDRO_615) -
              DENDRO_298 * (0.5 * DENDRO_417 + DENDRO_612) -
              DENDRO_298 * (DENDRO_426 + DENDRO_614) -
              DENDRO_298 * (DENDRO_165 * DENDRO_54 + DENDRO_183 + DENDRO_613) -
              DENDRO_606 * (DENDRO_168 + DENDRO_432 + DENDRO_433) -
              DENDRO_606 * (DENDRO_436 + DENDRO_437 + DENDRO_578) -
              DENDRO_630 * (DENDRO_327 + DENDRO_464 + DENDRO_520) +
              DENDRO_652)) -
    DENDRO_22 *
        (DENDRO_593 +
         alpha[pp] *
             (-DENDRO_104 * (DENDRO_551 + DENDRO_595) -
              DENDRO_104 * (DENDRO_552 + DENDRO_594) -
              DENDRO_104 * (DENDRO_210 * DENDRO_213 + DENDRO_389 + DENDRO_549) -
              DENDRO_114 * (0.5 * DENDRO_315 + DENDRO_547) -
              DENDRO_114 * (-DENDRO_337 + DENDRO_548) +
              DENDRO_135 * (DENDRO_563 + DENDRO_600) +
              DENDRO_135 * (DENDRO_564 + DENDRO_568) +
              DENDRO_135 * (DENDRO_569 - DENDRO_598) +
              DENDRO_135 * (-DENDRO_203 * DENDRO_213 + DENDRO_599) +
              DENDRO_135 * (DENDRO_210 * DENDRO_241 + DENDRO_562 + DENDRO_600) +
              DENDRO_156 * (DENDRO_376 + DENDRO_571) +
              DENDRO_156 * (DENDRO_379 + DENDRO_577) +
              DENDRO_156 * (DENDRO_379 + DENDRO_605) +
              DENDRO_156 * (DENDRO_443 + DENDRO_603) +
              DENDRO_156 * (DENDRO_444 + DENDRO_604) +
              DENDRO_156 * (DENDRO_574 + DENDRO_601) -
              DENDRO_180 * (DENDRO_533 + DENDRO_555) -
              DENDRO_180 * (DENDRO_558 + DENDRO_597) -
              DENDRO_180 * (DENDRO_561 + DENDRO_597) -
              DENDRO_298 * (DENDRO_416 + DENDRO_554) -
              DENDRO_298 * (DENDRO_171 * DENDRO_173 + DENDRO_596) +
              DENDRO_311 * (DENDRO_213 * grad_1_gt5[pp] + DENDRO_410) +
              DENDRO_588 - DENDRO_606 * (DENDRO_519 + DENDRO_556) +
              DENDRO_607)) -
    DENDRO_22 *
        (DENDRO_593 +
         alpha[pp] *
             (-DENDRO_104 * (-DENDRO_384 + DENDRO_551) -
              DENDRO_104 * (0.5 * DENDRO_388 + DENDRO_549) -
              DENDRO_104 * (DENDRO_200 * DENDRO_244 + DENDRO_386 + DENDRO_552) -
              DENDRO_114 * (1.0 * DENDRO_309 + DENDRO_546) -
              DENDRO_114 * (0.5 * DENDRO_330 + DENDRO_548) -
              DENDRO_114 * (DENDRO_200 * DENDRO_241 + DENDRO_321 + DENDRO_547) +
              DENDRO_135 * (DENDRO_566 + DENDRO_567) +
              DENDRO_135 * (DENDRO_568 + DENDRO_569) +
              DENDRO_135 * (-DENDRO_166 * DENDRO_224 + DENDRO_564) +
              DENDRO_135 * (DENDRO_200 * DENDRO_213 + DENDRO_566) +
              DENDRO_135 * (DENDRO_211 * DENDRO_241 + DENDRO_563) +
              DENDRO_140 * (DENDRO_487 + DENDRO_545) +
              DENDRO_140 * (DENDRO_543 + DENDRO_544) +
              DENDRO_156 * (DENDRO_444 + DENDRO_573) +
              DENDRO_156 * (DENDRO_570 + DENDRO_571) +
              DENDRO_156 * (DENDRO_573 + DENDRO_574) +
              DENDRO_156 * (DENDRO_175 * DENDRO_385 + DENDRO_577) -
              DENDRO_180 * (DENDRO_317 + DENDRO_467) -
              DENDRO_180 * (DENDRO_317 + DENDRO_559) -
              DENDRO_180 * (DENDRO_333 + DENDRO_555) -
              DENDRO_180 * (DENDRO_461 + DENDRO_557) -
              DENDRO_180 * (DENDRO_463 + DENDRO_558) -
              DENDRO_180 * (DENDRO_474 + DENDRO_561) -
              DENDRO_298 * (DENDRO_415 + DENDRO_554) -
              DENDRO_298 * (DENDRO_116 * DENDRO_165 + DENDRO_538) +
              DENDRO_311 * (DENDRO_241 * grad_2_gt3[pp] + DENDRO_292) -
              DENDRO_579 * (DENDRO_438 + DENDRO_578) + DENDRO_588)) +
    DENDRO_25 *
        (-DENDRO_266 * DENDRO_366 +
         DENDRO_270 * (DENDRO_227 + DENDRO_362 * DENDRO_72) +
         DENDRO_363 * (DENDRO_244 + DENDRO_273 * DENDRO_362) +
         alpha[pp] *
             (DENDRO_114 * DENDRO_237 * DENDRO_401 +
              DENDRO_135 * (DENDRO_386 - DENDRO_387) +
              DENDRO_135 * (DENDRO_389 + 1.0 * DENDRO_390) +
              DENDRO_135 * (DENDRO_227 * DENDRO_91 - DENDRO_384) +
              DENDRO_141 * (DENDRO_394 + DENDRO_396) +
              DENDRO_141 * (DENDRO_173 * grad_2_gt3[pp] + DENDRO_400) +
              DENDRO_141 * (DENDRO_222 * grad_0_gt3[pp] + DENDRO_397) +
              DENDRO_156 * (DENDRO_392 + DENDRO_393) +
              DENDRO_156 * (DENDRO_395 + 1.0 * DENDRO_396) +
              DENDRO_156 * (DENDRO_244 * DENDRO_91 + DENDRO_391) -
              DENDRO_173 * DENDRO_405 - DENDRO_180 * (DENDRO_374 + DENDRO_376) -
              DENDRO_180 * (-1.0 * DENDRO_380 - DENDRO_382) -
              DENDRO_180 * (DENDRO_222 * DENDRO_378 + DENDRO_377) -
              DENDRO_180 * (DENDRO_237 * DENDRO_378 + DENDRO_379) -
              DENDRO_191 *
                  (DENDRO_187 * DENDRO_369 +
                   DENDRO_188 * (DENDRO_370 + DENDRO_371 + DENDRO_372) +
                   DENDRO_190 * DENDRO_373 + DENDRO_367) -
              DENDRO_214 * DENDRO_294 * grad_1_gt3[pp] -
              DENDRO_221 * DENDRO_403 - DENDRO_222 * DENDRO_404 +
              DENDRO_311 * (DENDRO_388 + DENDRO_390) +
              DENDRO_311 * (DENDRO_221 * grad_0_gt3[pp] + DENDRO_398) +
              DENDRO_311 * (DENDRO_237 * grad_2_gt3[pp] + DENDRO_399) +
              DENDRO_355 * grad_1_gt3[pp] + DENDRO_357 * grad_0_gt3[pp] +
              DENDRO_359 * grad_2_gt3[pp] + DENDRO_361 * DENDRO_365 +
              DENDRO_412) -
         4 * grad2_1_1_alpha[pp]) +
    DENDRO_32 *
        (DENDRO_517 +
         alpha[pp] *
             (-DENDRO_104 * (DENDRO_407 + DENDRO_446) -
              DENDRO_104 * (DENDRO_221 * DENDRO_47 + DENDRO_443) -
              DENDRO_114 * (1.0 * DENDRO_312 + DENDRO_439) -
              DENDRO_114 * (0.5 * DENDRO_324 + DENDRO_441) -
              DENDRO_114 * (DENDRO_171 * DENDRO_241 + DENDRO_320 + DENDRO_440) +
              DENDRO_135 * (DENDRO_322 + DENDRO_461) +
              DENDRO_135 * (DENDRO_322 + DENDRO_469) +
              DENDRO_135 * (DENDRO_462 + DENDRO_463) +
              DENDRO_135 * (DENDRO_465 + DENDRO_467) +
              DENDRO_135 * (DENDRO_471 + DENDRO_474) +
              DENDRO_135 * (0.5 * DENDRO_326 + DENDRO_328 - DENDRO_459) +
              DENDRO_140 * (DENDRO_432 + DENDRO_435) +
              DENDRO_140 * (DENDRO_436 + DENDRO_438) +
              DENDRO_156 * (DENDRO_415 + DENDRO_480) +
              DENDRO_156 * (DENDRO_475 + DENDRO_476) +
              DENDRO_156 * (DENDRO_478 + DENDRO_481) +
              DENDRO_156 * (DENDRO_182 * DENDRO_237 + DENDRO_482 + DENDRO_484) -
              DENDRO_161 * (DENDRO_241 * grad_2_gt0[pp] + DENDRO_289) -
              DENDRO_180 * (DENDRO_454 + DENDRO_456) -
              DENDRO_180 * (DENDRO_457 + DENDRO_458) -
              DENDRO_180 * (DENDRO_204 * DENDRO_48 + DENDRO_451) -
              DENDRO_180 * (DENDRO_171 * DENDRO_186 + DENDRO_454 + DENDRO_455) -
              DENDRO_180 * (DENDRO_241 * DENDRO_51 + DENDRO_452 + DENDRO_453) -
              DENDRO_298 * (DENDRO_179 + DENDRO_448) -
              DENDRO_298 * (0.5 * DENDRO_419 + DENDRO_447) -
              DENDRO_298 * (DENDRO_171 * DENDRO_189 + DENDRO_427 + DENDRO_449) -
              DENDRO_488 * (DENDRO_485 + DENDRO_487) + DENDRO_509)) +
    DENDRO_32 *
        (DENDRO_517 +
         alpha[pp] *
             (-DENDRO_104 * (DENDRO_173 * DENDRO_200 + DENDRO_524) -
              DENDRO_114 * (0.5 * DENDRO_313 + DENDRO_440) -
              DENDRO_114 * (DENDRO_224 * DENDRO_50 + DENDRO_340 + DENDRO_439) +
              DENDRO_134 * DENDRO_522 + DENDRO_134 * (DENDRO_332 + DENDRO_519) -
              DENDRO_135 * DENDRO_531 + DENDRO_135 * (DENDRO_462 + DENDRO_532) +
              DENDRO_135 * (DENDRO_469 + DENDRO_533) +
              DENDRO_135 * (DENDRO_471 + DENDRO_532) + DENDRO_156 * DENDRO_534 +
              DENDRO_156 * DENDRO_541 + DENDRO_156 * (DENDRO_416 + DENDRO_481) +
              DENDRO_156 * (DENDRO_480 + DENDRO_535) +
              DENDRO_156 * (DENDRO_484 + DENDRO_537) +
              DENDRO_156 * (DENDRO_536 + DENDRO_537) -
              DENDRO_161 * (DENDRO_186 * grad_0_gt5[pp] + DENDRO_414) -
              DENDRO_180 * (DENDRO_453 + DENDRO_528) -
              DENDRO_180 * (-DENDRO_172 * DENDRO_186 + DENDRO_456) -
              DENDRO_180 * (DENDRO_241 * DENDRO_50 + DENDRO_528) -
              DENDRO_298 * (1.0 * DENDRO_421 + DENDRO_449) -
              DENDRO_298 * (DENDRO_186 * DENDRO_50 + DENDRO_423 + DENDRO_447) -
              DENDRO_488 * (DENDRO_375 + DENDRO_543) + DENDRO_509 - DENDRO_523 -
              DENDRO_525 - DENDRO_526 - DENDRO_527 - DENDRO_529)) +
    DENDRO_34 *
        (DENDRO_267 * DENDRO_55 - DENDRO_269 * DENDRO_73 -
         DENDRO_363 * DENDRO_61 +
         alpha[pp] *
             (-DENDRO_104 * DENDRO_124 * DENDRO_173 - DENDRO_105 +
              DENDRO_112 * DENDRO_114 * DENDRO_175 - DENDRO_121 -
              DENDRO_129 * DENDRO_131 * DENDRO_246 - DENDRO_130 +
              DENDRO_135 * DENDRO_152 +
              DENDRO_135 * (-1.0 * DENDRO_167 - DENDRO_170) +
              DENDRO_135 * (-1.0 * DENDRO_174 - DENDRO_177) +
              DENDRO_135 * (DENDRO_151 * DENDRO_175 + DENDRO_424) +
              DENDRO_135 * (DENDRO_218 * DENDRO_47 + DENDRO_416) +
              DENDRO_135 * (DENDRO_222 * DENDRO_50 + DENDRO_415) +
              DENDRO_141 * DENDRO_146 + DENDRO_141 * (DENDRO_417 + DENDRO_418) +
              DENDRO_141 * (DENDRO_173 * grad_2_gt0[pp] + DENDRO_425) +
              DENDRO_156 * DENDRO_184 +
              DENDRO_156 * (1.0 * DENDRO_418 + DENDRO_422) +
              DENDRO_156 * (DENDRO_150 * DENDRO_189 + DENDRO_426) -
              DENDRO_161 * (DENDRO_419 + DENDRO_420) -
              DENDRO_161 * (DENDRO_175 * grad_2_gt0[pp] + DENDRO_421) -
              DENDRO_162 - DENDRO_180 * (1.0 * DENDRO_420 + DENDRO_423) -
              DENDRO_180 * (-2 * DENDRO_172 * DENDRO_189 + DENDRO_427) -
              DENDRO_181 -
              DENDRO_191 * (DENDRO_185 +
                            DENDRO_187 * (DENDRO_428 + DENDRO_429 + DENDRO_66) +
                            DENDRO_188 * DENDRO_430 + DENDRO_190 * DENDRO_431) +
              DENDRO_217 * DENDRO_354 -
              6.0 * DENDRO_229 * DENDRO_83 * grad_0_gt0[pp] +
              DENDRO_232 * DENDRO_356 + DENDRO_249 * DENDRO_358 +
              DENDRO_262 * DENDRO_360 + DENDRO_276 * grad_0_Gt2[pp] +
              DENDRO_282 * grad2_0_2_gt0[pp] + DENDRO_284 * grad2_0_0_gt0[pp] +
              DENDRO_286 * grad2_1_1_gt0[pp] + DENDRO_288 * grad2_2_2_gt0[pp] -
              DENDRO_293 * DENDRO_413 + DENDRO_406 * grad_0_Gt1[pp] -
              DENDRO_411 * DENDRO_414 - DENDRO_77 - DENDRO_80 - DENDRO_82 +
              4 * grad_0_Gt0[pp] * gt0[pp]) -
         4 * grad2_0_0_alpha[pp]) +
    DENDRO_38 *
        (DENDRO_265 * DENDRO_267 - DENDRO_268 * DENDRO_270 -
         DENDRO_274 * DENDRO_275 +
         alpha[pp] *
             (-DENDRO_104 * DENDRO_221 * DENDRO_302 -
              DENDRO_129 * DENDRO_225 * DENDRO_303 + DENDRO_135 * DENDRO_339 +
              DENDRO_135 * (1.0 * DENDRO_316 + DENDRO_321) +
              DENDRO_135 * (DENDRO_224 * DENDRO_86 - DENDRO_337) +
              DENDRO_156 * DENDRO_319 + DENDRO_156 * DENDRO_329 +
              DENDRO_156 * (DENDRO_331 + DENDRO_333) +
              DENDRO_156 * (DENDRO_171 * DENDRO_237 + DENDRO_307) +
              DENDRO_156 * (DENDRO_175 * DENDRO_200 + DENDRO_305) +
              DENDRO_156 * (DENDRO_218 * DENDRO_318 + DENDRO_322) -
              DENDRO_161 * (DENDRO_313 + DENDRO_314) -
              DENDRO_161 * (DENDRO_218 * grad_0_gt5[pp] + DENDRO_312) -
              DENDRO_180 * (1.0 * DENDRO_314 + DENDRO_320) -
              DENDRO_180 * (DENDRO_204 * DENDRO_86 - DENDRO_335) -
              DENDRO_180 * (DENDRO_224 * DENDRO_341 + DENDRO_340) -
              DENDRO_191 *
                  (DENDRO_187 * DENDRO_347 + DENDRO_188 * DENDRO_349 +
                   DENDRO_190 * (DENDRO_351 + DENDRO_352 + DENDRO_353) +
                   DENDRO_342) -
              DENDRO_218 * DENDRO_297 * DENDRO_298 -
              DENDRO_242 * DENDRO_294 * grad_2_gt5[pp] +
              DENDRO_272 * DENDRO_361 + DENDRO_276 * grad_2_Gt0[pp] +
              DENDRO_277 * grad_2_Gt1[pp] - DENDRO_278 - DENDRO_279 -
              DENDRO_280 + DENDRO_282 * grad2_0_2_gt5[pp] +
              DENDRO_284 * grad2_0_0_gt5[pp] + DENDRO_286 * grad2_1_1_gt5[pp] +
              DENDRO_288 * grad2_2_2_gt5[pp] - DENDRO_289 * DENDRO_291 -
              DENDRO_292 * DENDRO_293 - DENDRO_296 - DENDRO_301 - DENDRO_304 +
              DENDRO_310 * DENDRO_311 + DENDRO_311 * (DENDRO_315 + DENDRO_316) +
              DENDRO_311 * (DENDRO_221 * grad_0_gt5[pp] + DENDRO_330) -
              DENDRO_325 + DENDRO_355 * grad_1_gt5[pp] +
              DENDRO_357 * grad_0_gt5[pp] + DENDRO_359 * grad_2_gt5[pp] +
              4 * grad_2_Gt2[pp] * gt5[pp]) -
         4 * grad2_2_2_alpha[pp]);
const double DENDRO_666 = DENDRO_30 * DENDRO_665;
const double DENDRO_667 = (1.0 / 12.0) * chi[pp];
const double DENDRO_668 = (1.0 / 3.0) * At1[pp];
const double DENDRO_669 = At1[pp] * DENDRO_20;
const double DENDRO_670 = At4[pp] * DENDRO_22;
const double DENDRO_671 = At3[pp] * DENDRO_25;
const double DENDRO_672 = DENDRO_669 + DENDRO_670 - DENDRO_671;
const double DENDRO_673 = At4[pp] * DENDRO_32;
const double DENDRO_674 =
    -At1[pp] * DENDRO_34 + At3[pp] * DENDRO_20 - DENDRO_673;
const double DENDRO_675 = At4[pp] * DENDRO_38;
const double DENDRO_676 =
    -At1[pp] * DENDRO_32 + At3[pp] * DENDRO_22 - DENDRO_675;
const double DENDRO_677 = 6.0 * grad_2_alpha[pp];
const double DENDRO_678 = 6.0 * grad_0_alpha[pp];
const double DENDRO_679 = 6.0 * grad_1_alpha[pp];
const double DENDRO_680 = DENDRO_91 * DENDRO_93;
const double DENDRO_681 = DENDRO_220 * grad_2_gt0[pp];
const double DENDRO_682 = DENDRO_150 * DENDRO_93;
const double DENDRO_683 = DENDRO_681 + DENDRO_682;
const double DENDRO_684 = DENDRO_383 * DENDRO_88;
const double DENDRO_685 = DENDRO_127 * DENDRO_385;
const double DENDRO_686 = DENDRO_178 * DENDRO_93;
const double DENDRO_687 =
    DENDRO_151 * DENDRO_94 + 0.25 * DENDRO_220 * grad_0_gt0[pp];
const double DENDRO_688 = DENDRO_133 + DENDRO_137;
const double DENDRO_689 = DENDRO_108 * DENDRO_575;
const double DENDRO_690 = DENDRO_200 * DENDRO_60;
const double DENDRO_691 = DENDRO_136 * DENDRO_236 + DENDRO_690;
const double DENDRO_692 = DENDRO_150 * DENDRO_88;
const double DENDRO_693 = DENDRO_220 * grad_1_gt0[pp] + DENDRO_692;
const double DENDRO_694 = 1.0 * DENDRO_134;
const double DENDRO_695 = DENDRO_150 * DENDRO_236;
const double DENDRO_696 = 2.0 * DENDRO_231;
const double DENDRO_697 = 2.0 * DENDRO_216;
const double DENDRO_698 = 2.0 * DENDRO_248;
const double DENDRO_699 = DENDRO_261 * DENDRO_44;
const double DENDRO_700 = (1.0 / 3.0) * At2[pp];
const double DENDRO_701 = At5[pp] * DENDRO_22;
const double DENDRO_702 =
    At2[pp] * DENDRO_20 - At4[pp] * DENDRO_25 + DENDRO_701;
const double DENDRO_703 = At5[pp] * DENDRO_32;
const double DENDRO_704 =
    -At2[pp] * DENDRO_34 + At4[pp] * DENDRO_20 - DENDRO_703;
const double DENDRO_705 = At4[pp] * DENDRO_22 - At5[pp] * DENDRO_38 - DENDRO_33;
const double DENDRO_706 = DENDRO_108 * grad_2_gt5[pp];
const double DENDRO_707 = DENDRO_88 * grad_0_gt5[pp];
const double DENDRO_708 = DENDRO_108 * DENDRO_306;
const double DENDRO_709 = 0.25 * DENDRO_127 * grad_2_gt5[pp];
const double DENDRO_710 = DENDRO_708 + DENDRO_709;
const double DENDRO_711 = DENDRO_86 * DENDRO_88;
const double DENDRO_712 = DENDRO_108 * DENDRO_109 + DENDRO_198 * DENDRO_60;
const double DENDRO_713 = -0.5 * DENDRO_172 * DENDRO_93;
const double DENDRO_714 = -0.5 * DENDRO_172 * DENDRO_88;
const double DENDRO_715 = DENDRO_182 * DENDRO_88;
const double DENDRO_716 = 2 * At4[pp];
const double DENDRO_717 = At3[pp] * DENDRO_36;
const double DENDRO_718 = DENDRO_30 * DENDRO_716;
const double DENDRO_719 = 0.5 * DENDRO_362;
const double DENDRO_720 = DENDRO_30 * DENDRO_74;
const double DENDRO_721 = DENDRO_212 * grad_2_gt3[pp];
const double DENDRO_722 = DENDRO_127 * DENDRO_306;
const double DENDRO_723 = 1.0 * DENDRO_220;
const double DENDRO_724 = 0.25 * DENDRO_692;
const double DENDRO_725 = DENDRO_212 * grad_0_gt3[pp];
const double DENDRO_726 = (1.0 / 3.0) * At4[pp];
const double DENDRO_727 = DENDRO_236 * grad_2_gt5[pp];
const double DENDRO_728 = DENDRO_109 * DENDRO_236;
const double DENDRO_729 = DENDRO_709 + DENDRO_728;
const double DENDRO_730 = DENDRO_236 * DENDRO_306;
const double DENDRO_731 = -0.5 * DENDRO_244 * grad_0_gt5[pp] + DENDRO_722;
const double DENDRO_732 = DENDRO_240 * grad_0_gt5[pp];
const double DENDRO_733 = DENDRO_240 * grad_1_gt5[pp];
const double DENDRO_734 = DENDRO_20 * grad_1_chi[pp] +
                          DENDRO_343 * grad_2_chi[pp] +
                          DENDRO_345 * grad_0_chi[pp];
const double DENDRO_735 = 0.5 * DENDRO_264;
const double DENDRO_736 = DENDRO_30 * grad_0_alpha[pp];
const double DENDRO_737 = DENDRO_348 * grad_1_chi[pp] + DENDRO_41;
const double DENDRO_738 = DENDRO_30 * grad_1_alpha[pp];
const double DENDRO_739 = DENDRO_22 * grad_1_chi[pp] +
                          DENDRO_343 * grad_0_chi[pp] +
                          DENDRO_350 * grad_2_chi[pp];
const double DENDRO_740 = 0.5 * DENDRO_739;
const double DENDRO_741 = DENDRO_30 * grad_2_alpha[pp];
const double DENDRO_742 = 0.5 * DENDRO_737;
const double DENDRO_743 = 0.5 * grad_1_alpha[pp];
const double DENDRO_744 = 0.5 * grad_2_alpha[pp];
const double DENDRO_745 = 0.5 * grad_0_alpha[pp];
const double DENDRO_746 = DENDRO_254 * DENDRO_30;
const double DENDRO_747 = (DENDRO_20 * DENDRO_20);
const double DENDRO_748 = (DENDRO_32 * DENDRO_32);
const double DENDRO_749 = 2 * DENDRO_34;
const double DENDRO_750 = At0[pp] * (DENDRO_34 * DENDRO_34) +
                          At3[pp] * DENDRO_747 + At5[pp] * DENDRO_748 -
                          DENDRO_252 * DENDRO_673 + DENDRO_33 * DENDRO_749 -
                          DENDRO_669 * DENDRO_749;
const double DENDRO_751 = 3 * DENDRO_83;
const double DENDRO_752 = (DENDRO_22 * DENDRO_22);
const double DENDRO_753 = 2 * DENDRO_25;
const double DENDRO_754 = At0[pp] * DENDRO_747 +
                          At3[pp] * (DENDRO_25 * DENDRO_25) +
                          At5[pp] * DENDRO_752 + DENDRO_23 * DENDRO_252 -
                          DENDRO_669 * DENDRO_753 - DENDRO_670 * DENDRO_753;
const double DENDRO_755 = At1[pp] * DENDRO_22;
const double DENDRO_756 = 2 * DENDRO_38;
const double DENDRO_757 = At0[pp] * DENDRO_748 + At3[pp] * DENDRO_752 +
                          At5[pp] * (DENDRO_38 * DENDRO_38) -
                          DENDRO_254 * DENDRO_755 + DENDRO_33 * DENDRO_756 -
                          DENDRO_670 * DENDRO_756;
const double DENDRO_758 = At3[pp] * DENDRO_20;
const double DENDRO_759 =
    At2[pp] * DENDRO_748 - DENDRO_20 * DENDRO_675 - DENDRO_22 * DENDRO_673 +
    DENDRO_22 * DENDRO_758 + DENDRO_32 * DENDRO_35 - DENDRO_32 * DENDRO_669 +
    DENDRO_34 * DENDRO_39 - DENDRO_34 * DENDRO_755 + DENDRO_38 * DENDRO_703;
const double DENDRO_760 = 6 * DENDRO_83;
const double DENDRO_761 =
    -At1[pp] * DENDRO_25 * DENDRO_34 - At1[pp] * DENDRO_747 -
    At4[pp] * DENDRO_20 * DENDRO_22 - At4[pp] * DENDRO_25 * DENDRO_32 +
    DENDRO_20 * DENDRO_33 + DENDRO_20 * DENDRO_35 + DENDRO_22 * DENDRO_703 +
    DENDRO_23 * DENDRO_34 + DENDRO_25 * DENDRO_758;
const double DENDRO_762 = -DENDRO_761;
const double DENDRO_763 =
    -At1[pp] * DENDRO_20 * DENDRO_22 - At1[pp] * DENDRO_25 * DENDRO_32 -
    At4[pp] * DENDRO_25 * DENDRO_38 - At4[pp] * DENDRO_752 +
    DENDRO_20 * DENDRO_39 + DENDRO_21 * DENDRO_32 + DENDRO_22 * DENDRO_33 +
    DENDRO_22 * DENDRO_671 + DENDRO_38 * DENDRO_701;
const double DENDRO_764 = -DENDRO_763;
const double DENDRO_765 = (1.0 / 3.0) * alpha[pp];
const double DENDRO_766 = (7.0 / 3.0) * DENDRO_281;
const double DENDRO_767 = (1.0 / 3.0) * DENDRO_281;
const double DENDRO_768 = (1.0 / 3.0) * DENDRO_283;
const double DENDRO_769 = 2 * DENDRO_83;
const double DENDRO_770 = DENDRO_769 * grad_0_alpha[pp];
const double DENDRO_771 = pow(DENDRO_29, -3);
const double DENDRO_772 = 4 * grad_0_K[pp];
const double DENDRO_773 = DENDRO_30 * DENDRO_765;
const double DENDRO_774 = DENDRO_769 * grad_2_alpha[pp];
const double DENDRO_775 = DENDRO_769 * grad_1_alpha[pp];
const double DENDRO_776 = 4 * grad_2_K[pp];
const double DENDRO_777 = 4 * grad_1_K[pp];
const double DENDRO_778 = 9 * DENDRO_44;
const double DENDRO_779 = DENDRO_188 * DENDRO_778;
const double DENDRO_780 =
    -2.0 / 3.0 * DENDRO_16 * DENDRO_259 * DENDRO_83 -
    2 * DENDRO_186 * DENDRO_750 * DENDRO_771 * alpha[pp] -
    7.0 / 3.0 * DENDRO_20 * DENDRO_30 * grad2_0_1_beta0[pp] -
    1.0 / 3.0 * DENDRO_20 * DENDRO_30 * grad2_1_1_beta1[pp] -
    1.0 / 3.0 * DENDRO_20 * DENDRO_30 * grad2_1_2_beta2[pp] -
    2.0 * DENDRO_218 * DENDRO_759 * DENDRO_771 * alpha[pp] -
    2 * DENDRO_22 * DENDRO_30 * grad2_1_2_beta0[pp] -
    2.0 * DENDRO_221 * DENDRO_764 * DENDRO_771 * alpha[pp] -
    2.0 * DENDRO_222 * DENDRO_762 * DENDRO_771 * alpha[pp] -
    2 * DENDRO_224 * DENDRO_757 * DENDRO_771 * alpha[pp] -
    2 * DENDRO_227 * DENDRO_754 * DENDRO_771 * alpha[pp] +
    (4.0 / 3.0) * DENDRO_283 * grad2_0_0_beta0[pp] +
    DENDRO_285 * grad2_1_1_beta0[pp] + DENDRO_287 * grad2_2_2_beta0[pp] +
    DENDRO_354 * grad_1_beta0[pp] + DENDRO_356 * grad_0_beta0[pp] +
    DENDRO_358 * grad_2_beta0[pp] + DENDRO_750 * DENDRO_770 +
    DENDRO_759 * DENDRO_774 + DENDRO_762 * DENDRO_775 +
    DENDRO_766 * grad2_0_2_beta0[pp] + DENDRO_767 * grad2_1_2_beta1[pp] +
    DENDRO_767 * grad2_2_2_beta2[pp] + DENDRO_768 * grad2_0_1_beta1[pp] +
    DENDRO_768 * grad2_0_2_beta2[pp] +
    DENDRO_773 * (DENDRO_20 * DENDRO_777 + DENDRO_762 * DENDRO_779) +
    DENDRO_773 * (9 * DENDRO_30 * DENDRO_44 * DENDRO_750 * grad_0_chi[pp] -
                  DENDRO_34 * DENDRO_772) +
    DENDRO_773 * (9 * DENDRO_30 * DENDRO_44 * DENDRO_759 * grad_2_chi[pp] -
                  DENDRO_32 * DENDRO_776) -
    beta0[pp] * grad_0_Gt0[pp] - beta1[pp] * grad_1_Gt0[pp] -
    beta2[pp] * grad_2_Gt0[pp];
const double DENDRO_781 = DENDRO_283 * grad2_0_0_beta1[pp];
const double DENDRO_782 = DENDRO_287 * grad2_2_2_beta1[pp];
const double DENDRO_783 = DENDRO_746 * grad2_0_2_beta1[pp];
const double DENDRO_784 = DENDRO_754 * DENDRO_775;
const double DENDRO_785 = (4.0 / 3.0) * DENDRO_285 * grad2_1_1_beta1[pp];
const double DENDRO_786 = DENDRO_20 * DENDRO_772;
const double DENDRO_787 = DENDRO_187 * DENDRO_778;
const double DENDRO_788 = DENDRO_22 * DENDRO_776;
const double DENDRO_789 = DENDRO_190 * DENDRO_778;
const double DENDRO_790 = (1.0 / 3.0) * DENDRO_285;
const double DENDRO_791 = DENDRO_790 * grad2_0_1_beta0[pp];
const double DENDRO_792 = DENDRO_790 * grad2_1_2_beta2[pp];
const double DENDRO_793 = DENDRO_25 * DENDRO_777 - DENDRO_754 * DENDRO_779;
const double DENDRO_794 = DENDRO_20 * DENDRO_30;
const double DENDRO_795 = (1.0 / 3.0) * DENDRO_794;
const double DENDRO_796 = DENDRO_22 * DENDRO_30;
const double DENDRO_797 = (1.0 / 3.0) * grad2_0_2_beta0[pp];
const double DENDRO_798 = DENDRO_1 * DENDRO_771;
const double DENDRO_799 = 2.0 * DENDRO_771 * alpha[pp];
const double DENDRO_800 = beta0[pp] * grad_0_Gt1[pp] +
                          beta1[pp] * grad_1_Gt1[pp] +
                          beta2[pp] * grad_2_Gt1[pp];
const double DENDRO_801 =
    -2.0 / 3.0 * DENDRO_16 * DENDRO_260 * DENDRO_83 -
    2.0 * DENDRO_173 * DENDRO_762 * DENDRO_771 * alpha[pp] -
    2.0 * DENDRO_175 * DENDRO_759 * DENDRO_771 * alpha[pp] -
    2 * DENDRO_189 * DENDRO_750 * DENDRO_771 * alpha[pp] -
    2 * DENDRO_20 * DENDRO_30 * grad2_0_1_beta2[pp] -
    1.0 / 3.0 * DENDRO_22 * DENDRO_30 * grad2_0_1_beta0[pp] -
    1.0 / 3.0 * DENDRO_22 * DENDRO_30 * grad2_1_1_beta1[pp] -
    7.0 / 3.0 * DENDRO_22 * DENDRO_30 * grad2_1_2_beta2[pp] -
    2.0 * DENDRO_237 * DENDRO_764 * DENDRO_771 * alpha[pp] -
    2 * DENDRO_241 * DENDRO_757 * DENDRO_771 * alpha[pp] -
    2 * DENDRO_244 * DENDRO_754 * DENDRO_771 * alpha[pp] +
    DENDRO_283 * grad2_0_0_beta2[pp] + DENDRO_285 * grad2_1_1_beta2[pp] +
    DENDRO_287 * DENDRO_797 + (1.0 / 3.0) * DENDRO_287 * grad2_1_2_beta1[pp] +
    (4.0 / 3.0) * DENDRO_287 * grad2_2_2_beta2[pp] +
    DENDRO_354 * grad_1_beta2[pp] + DENDRO_356 * grad_0_beta2[pp] +
    DENDRO_358 * grad_2_beta2[pp] + DENDRO_757 * DENDRO_774 +
    DENDRO_759 * DENDRO_770 + DENDRO_764 * DENDRO_775 +
    DENDRO_766 * grad2_0_2_beta2[pp] + DENDRO_767 * grad2_0_0_beta0[pp] +
    DENDRO_767 * grad2_0_1_beta1[pp] +
    DENDRO_773 * (DENDRO_22 * DENDRO_777 + DENDRO_764 * DENDRO_779) +
    DENDRO_773 * (9 * DENDRO_30 * DENDRO_44 * DENDRO_757 * grad_2_chi[pp] -
                  DENDRO_38 * DENDRO_776) +
    DENDRO_773 * (9 * DENDRO_30 * DENDRO_44 * DENDRO_759 * grad_0_chi[pp] -
                  DENDRO_32 * DENDRO_772) -
    beta0[pp] * grad_0_Gt2[pp] - beta1[pp] * grad_1_Gt2[pp] -
    beta2[pp] * grad_2_Gt2[pp];

// Dendro: printing variables
//--
a_rhs[pp] =
    -K[pp] * (A_lambda[0] * pow(alpha[pp], 2) + A_lambda[1] * alpha[pp] +
              A_lambda[2]) +
    lambda[0] * (beta0[pp] * grad_0_alpha[pp] + beta1[pp] * grad_1_alpha[pp] +
                 beta2[pp] * grad_2_alpha[pp]);
//--
b_rhs0[pp]   = B0[pp] * DENDRO_0 + lambda[1] * (beta0[pp] * grad_0_beta0[pp] +
                                              beta1[pp] * grad_1_beta0[pp] +
                                              beta2[pp] * grad_2_beta0[pp]);
//--
b_rhs1[pp]   = B1[pp] * DENDRO_0 + lambda[1] * (beta0[pp] * grad_0_beta1[pp] +
                                              beta1[pp] * grad_1_beta1[pp] +
                                              beta2[pp] * grad_2_beta1[pp]);
//--
b_rhs2[pp]   = B2[pp] * DENDRO_0 + lambda[1] * (beta0[pp] * grad_0_beta2[pp] +
                                              beta1[pp] * grad_1_beta2[pp] +
                                              beta2[pp] * grad_2_beta2[pp]);
//--
gt_rhs00[pp] = -At0[pp] * DENDRO_1 + DENDRO_2 * gt0[pp] +
               DENDRO_3 * grad_0_beta1[pp] + DENDRO_4 * grad_0_beta2[pp] -
               DENDRO_5 * grad_1_beta1[pp] - DENDRO_5 * grad_2_beta2[pp] +
               beta0[pp] * grad_0_gt0[pp] + beta1[pp] * grad_1_gt0[pp] +
               beta2[pp] * grad_2_gt0[pp];
//--
gt_rhs01[pp] = -At1[pp] * DENDRO_1 + DENDRO_6 * grad_0_beta0[pp] +
               DENDRO_6 * grad_1_beta1[pp] - DENDRO_7 * gt1[pp] +
               beta0[pp] * grad_0_gt1[pp] + beta1[pp] * grad_1_gt1[pp] +
               beta2[pp] * grad_2_gt1[pp] + grad_0_beta1[pp] * gt3[pp] +
               grad_0_beta2[pp] * gt4[pp] + grad_1_beta0[pp] * gt0[pp] +
               grad_1_beta2[pp] * gt2[pp];
//--
gt_rhs02[pp] = -At2[pp] * DENDRO_1 + DENDRO_8 * grad_0_beta0[pp] +
               DENDRO_8 * grad_2_beta2[pp] - DENDRO_9 * gt2[pp] +
               beta0[pp] * grad_0_gt2[pp] + beta1[pp] * grad_1_gt2[pp] +
               beta2[pp] * grad_2_gt2[pp] + grad_0_beta1[pp] * gt4[pp] +
               grad_0_beta2[pp] * gt5[pp] + grad_2_beta0[pp] * gt0[pp] +
               grad_2_beta1[pp] * gt1[pp];
//--
gt_rhs11[pp] = -At3[pp] * DENDRO_1 - DENDRO_10 * gt3[pp] + DENDRO_11 * gt3[pp] +
               DENDRO_12 * grad_1_beta2[pp] + DENDRO_3 * grad_1_beta0[pp] -
               DENDRO_7 * gt3[pp] + beta0[pp] * grad_0_gt3[pp] +
               beta1[pp] * grad_1_gt3[pp] + beta2[pp] * grad_2_gt3[pp];
//--
gt_rhs12[pp] = -At4[pp] * DENDRO_1 - DENDRO_10 * gt4[pp] +
               DENDRO_13 * grad_1_beta1[pp] + DENDRO_13 * grad_2_beta2[pp] +
               beta0[pp] * grad_0_gt4[pp] + beta1[pp] * grad_1_gt4[pp] +
               beta2[pp] * grad_2_gt4[pp] + grad_1_beta0[pp] * gt2[pp] +
               grad_1_beta2[pp] * gt5[pp] + grad_2_beta0[pp] * gt1[pp] +
               grad_2_beta1[pp] * gt3[pp];
//--
gt_rhs22[pp] = -At5[pp] * DENDRO_1 - DENDRO_10 * gt5[pp] +
               DENDRO_12 * grad_2_beta1[pp] + DENDRO_14 * gt5[pp] +
               DENDRO_4 * grad_2_beta0[pp] - DENDRO_9 * gt5[pp] +
               beta0[pp] * grad_0_gt5[pp] + beta1[pp] * grad_1_gt5[pp] +
               beta2[pp] * grad_2_gt5[pp];
//--
chi_rhs[pp] = -DENDRO_15 * DENDRO_16 + DENDRO_15 * K[pp] * alpha[pp] +
              beta0[pp] * grad_0_chi[pp] + beta1[pp] * grad_1_chi[pp] +
              beta2[pp] * grad_2_chi[pp];
//--
At_rhs00[pp] =
    At0[pp] * DENDRO_2 - At0[pp] * DENDRO_7 - At0[pp] * DENDRO_9 +
    DENDRO_17 * grad_0_beta1[pp] + DENDRO_18 * grad_0_beta2[pp] +
    DENDRO_667 *
        (DENDRO_263 *
             (-DENDRO_105 - DENDRO_108 * DENDRO_112 * DENDRO_114 - DENDRO_121 +
              4 * DENDRO_124 * DENDRO_127 * DENDRO_25 * DENDRO_83 - DENDRO_130 +
              4 * DENDRO_131 * DENDRO_34 * DENDRO_60 * DENDRO_83 -
              DENDRO_135 * (DENDRO_133 + DENDRO_50 * DENDRO_88) -
              DENDRO_135 * (DENDRO_137 + DENDRO_47 * DENDRO_93) -
              DENDRO_135 * (DENDRO_167 + DENDRO_170) -
              DENDRO_135 * (DENDRO_174 + DENDRO_177) -
              DENDRO_135 * (DENDRO_150 * DENDRO_154 + DENDRO_153) -
              DENDRO_141 * (DENDRO_138 + DENDRO_139) -
              DENDRO_141 * (DENDRO_163 + DENDRO_164) +
              2.0 * DENDRO_146 * DENDRO_20 * DENDRO_83 +
              4 * DENDRO_152 * DENDRO_22 * DENDRO_83 -
              DENDRO_156 * (1.0 * DENDRO_139 + DENDRO_155) -
              DENDRO_156 *
                  (DENDRO_127 * DENDRO_182 + 1.0 * DENDRO_150 * DENDRO_60) -
              DENDRO_162 -
              DENDRO_180 *
                  (-DENDRO_154 * DENDRO_51 + 2 * DENDRO_172 * DENDRO_60) -
              DENDRO_181 + 4 * DENDRO_184 * DENDRO_20 * DENDRO_83 -
              DENDRO_191 * (DENDRO_185 + DENDRO_186 * DENDRO_187 +
                            DENDRO_188 * DENDRO_54 + DENDRO_189 * DENDRO_190) -
              DENDRO_216 * DENDRO_217 - DENDRO_231 * DENDRO_232 -
              DENDRO_248 * DENDRO_249 +
              2.0 * DENDRO_25 * DENDRO_30 * grad2_1_1_gt0[pp] +
              3.0 * DENDRO_25 * DENDRO_83 * DENDRO_88 * grad_1_gt0[pp] -
              DENDRO_261 * DENDRO_262 +
              4 * DENDRO_30 * DENDRO_32 * grad2_0_2_gt0[pp] +
              2.0 * DENDRO_30 * DENDRO_34 * grad2_0_0_gt0[pp] +
              2.0 * DENDRO_30 * DENDRO_38 * grad2_2_2_gt0[pp] +
              2.0 * DENDRO_32 * DENDRO_83 * (DENDRO_142 + DENDRO_143) +
              4 * DENDRO_32 * DENDRO_83 * (1.0 * DENDRO_143 + DENDRO_157) +
              2.0 * DENDRO_32 * DENDRO_83 * (DENDRO_147 + DENDRO_148) +
              6.0 * DENDRO_34 * DENDRO_83 * DENDRO_94 * grad_0_gt0[pp] +
              3.0 * DENDRO_38 * DENDRO_83 * DENDRO_93 * grad_2_gt0[pp] -
              DENDRO_77 - DENDRO_80 - DENDRO_82 + 4 * grad_0_Gt0[pp] * gt0[pp] +
              4 * grad_0_Gt1[pp] * gt1[pp] + 4 * grad_0_Gt2[pp] * gt2[pp]) +
         DENDRO_55 * DENDRO_57 - DENDRO_61 * DENDRO_63 + DENDRO_666 * gt0[pp] -
         DENDRO_73 * DENDRO_74 - 12 * grad2_0_0_alpha[pp]) -
    alpha[pp] *
        (-At0[pp] * K[pp] +
         DENDRO_31 * (-At1[pp] * DENDRO_25 + DENDRO_21 + DENDRO_23) +
         DENDRO_37 * (At1[pp] * DENDRO_20 - DENDRO_33 - DENDRO_35) +
         DENDRO_40 * (-At0[pp] * DENDRO_32 + At1[pp] * DENDRO_22 - DENDRO_39)) +
    beta0[pp] * grad_0_At0[pp] + beta1[pp] * grad_1_At0[pp] +
    beta2[pp] * grad_2_At0[pp];
//--
At_rhs01[pp] =
    At0[pp] * grad_1_beta0[pp] - At1[pp] * DENDRO_7 +
    At2[pp] * grad_1_beta2[pp] + At3[pp] * grad_0_beta1[pp] +
    At4[pp] * grad_0_beta2[pp] +
    DENDRO_667 *
        (DENDRO_263 *
             (-DENDRO_104 * (DENDRO_166 * DENDRO_212 + DENDRO_658) -
              DENDRO_104 * (-DENDRO_100 * DENDRO_88 + DENDRO_659 + DENDRO_684) -
              DENDRO_104 * (-DENDRO_127 * DENDRO_295 +
                            0.5 * DENDRO_150 * DENDRO_244 - DENDRO_685) +
              DENDRO_113 * (DENDRO_680 + DENDRO_683) + DENDRO_114 * DENDRO_657 -
              DENDRO_135 * DENDRO_661 +
              DENDRO_135 * (-DENDRO_151 * DENDRO_212 + DENDRO_662) +
              DENDRO_135 * (DENDRO_383 * DENDRO_93 - DENDRO_575 * DENDRO_88 +
                            DENDRO_602) -
              DENDRO_141 *
                  (DENDRO_88 * grad_1_gt0[pp] + DENDRO_94 * grad_0_gt3[pp]) +
              DENDRO_156 * (-DENDRO_212 * DENDRO_47 + DENDRO_663) -
              DENDRO_156 * (DENDRO_127 * DENDRO_468 + DENDRO_127 * DENDRO_483 +
                            DENDRO_210 * DENDRO_60) +
              DENDRO_156 * (-DENDRO_127 * DENDRO_575 - DENDRO_211 * DENDRO_60 +
                            0.5 * DENDRO_244 * grad_2_gt0[pp]) +
              DENDRO_156 * (DENDRO_166 * DENDRO_94 - DENDRO_178 * DENDRO_88 +
                            DENDRO_626) +
              DENDRO_180 * DENDRO_660 + DENDRO_180 * (DENDRO_686 + DENDRO_687) +
              DENDRO_180 * (DENDRO_689 + DENDRO_691) +
              DENDRO_180 * (DENDRO_378 * DENDRO_94 + DENDRO_688) -
              DENDRO_191 * (DENDRO_173 * DENDRO_498 + DENDRO_222 * DENDRO_492 +
                            DENDRO_494 * DENDRO_99 + DENDRO_631) +
              DENDRO_298 *
                  (DENDRO_155 + DENDRO_47 * DENDRO_94 + DENDRO_48 * DENDRO_94) +
              DENDRO_298 * (0.25 * DENDRO_163 + 0.5 * DENDRO_164 +
                            DENDRO_378 * DENDRO_60) +
              DENDRO_638 + DENDRO_639 + DENDRO_640 + DENDRO_641 + DENDRO_642 +
              DENDRO_643 + DENDRO_644 + DENDRO_645 + DENDRO_646 + DENDRO_647 +
              DENDRO_648 + DENDRO_649 + DENDRO_650 - DENDRO_651 * DENDRO_699 +
              DENDRO_664 -
              DENDRO_694 * (DENDRO_693 + DENDRO_93 * grad_0_gt3[pp]) -
              DENDRO_694 * (DENDRO_108 * grad_2_gt3[pp] +
                            DENDRO_127 * grad_1_gt5[pp] + DENDRO_695) -
              DENDRO_696 * grad_0_gt1[pp] - DENDRO_697 * grad_1_gt1[pp] -
              DENDRO_698 * grad_2_gt1[pp]) +
         DENDRO_30 * DENDRO_677 * (-DENDRO_127 - DENDRO_58 * DENDRO_655) +
         DENDRO_651 * DENDRO_665 - DENDRO_653 * DENDRO_678 -
         DENDRO_654 * DENDRO_679 - 12 * grad2_0_1_alpha[pp]) +
    DENDRO_668 * grad_0_beta0[pp] + DENDRO_668 * grad_1_beta1[pp] -
    alpha[pp] * (-At1[pp] * K[pp] + DENDRO_31 * DENDRO_672 +
                 DENDRO_37 * DENDRO_674 + DENDRO_40 * DENDRO_676) +
    beta0[pp] * grad_0_At1[pp] + beta1[pp] * grad_1_At1[pp] +
    beta2[pp] * grad_2_At1[pp];
//--
At_rhs02[pp] =
    At0[pp] * grad_2_beta0[pp] + At1[pp] * grad_2_beta1[pp] -
    At2[pp] * DENDRO_9 + At4[pp] * grad_0_beta1[pp] +
    At5[pp] * grad_0_beta2[pp] +
    DENDRO_667 *
        (DENDRO_263 *
             (-DENDRO_114 * (DENDRO_172 * DENDRO_240 - 0.5 * DENDRO_706) -
              DENDRO_114 * (-DENDRO_109 * DENDRO_93 - DENDRO_223 * DENDRO_50 -
                            DENDRO_713) -
              DENDRO_135 * DENDRO_531 -
              DENDRO_135 * (DENDRO_151 * DENDRO_240 + DENDRO_710) -
              DENDRO_156 * (DENDRO_153 + DENDRO_691) -
              DENDRO_156 * (DENDRO_687 + DENDRO_715) -
              DENDRO_156 * (DENDRO_318 * DENDRO_94 + DENDRO_688) -
              DENDRO_156 * (DENDRO_108 * DENDRO_483 + DENDRO_153 + DENDRO_690) -
              DENDRO_180 * (0.5 * DENDRO_108 * DENDRO_172 - DENDRO_712) -
              DENDRO_180 * (DENDRO_172 * DENDRO_94 - DENDRO_182 * DENDRO_93 -
                            DENDRO_223 * DENDRO_52) -
              DENDRO_191 * (DENDRO_116 * DENDRO_494 + DENDRO_175 * DENDRO_498 +
                            DENDRO_218 * DENDRO_492 + DENDRO_489) +
              4 * DENDRO_20 * DENDRO_534 * DENDRO_83 +
              4 * DENDRO_20 * DENDRO_541 * DENDRO_83 +
              DENDRO_22 * DENDRO_522 * DENDRO_83 +
              4 * DENDRO_22 * DENDRO_83 *
                  (0.5 * DENDRO_172 * DENDRO_236 - DENDRO_710) +
              4 * DENDRO_22 * DENDRO_83 *
                  (-DENDRO_223 * DENDRO_47 - DENDRO_468 * DENDRO_93 -
                   DENDRO_714) +
              2.0 * DENDRO_25 * DENDRO_30 * grad2_1_1_gt2[pp] +
              DENDRO_25 * DENDRO_83 * (DENDRO_693 + DENDRO_711) +
              4 * DENDRO_25 * DENDRO_83 *
                  (DENDRO_127 * DENDRO_200 + 0.25 * DENDRO_695) +
              4 * DENDRO_30 * DENDRO_32 * grad2_0_2_gt2[pp] +
              2.0 * DENDRO_30 * DENDRO_34 * grad2_0_0_gt2[pp] +
              2.0 * DENDRO_30 * DENDRO_38 * grad2_2_2_gt2[pp] +
              4 * DENDRO_32 * DENDRO_83 *
                  (DENDRO_240 * DENDRO_50 + DENDRO_712) +
              2.0 * DENDRO_32 * DENDRO_83 *
                  (DENDRO_93 * grad_2_gt0[pp] + DENDRO_94 * grad_0_gt5[pp]) +
              4 * DENDRO_34 * DENDRO_83 *
                  (0.25 * DENDRO_147 + 1.0 * DENDRO_148) +
              4 * DENDRO_34 * DENDRO_83 *
                  (DENDRO_157 + DENDRO_50 * DENDRO_94 + DENDRO_51 * DENDRO_94) -
              DENDRO_505 - DENDRO_506 - DENDRO_507 - DENDRO_508 * DENDRO_699 -
              DENDRO_523 - DENDRO_525 - DENDRO_526 - DENDRO_527 - DENDRO_529 -
              DENDRO_694 * (DENDRO_683 + DENDRO_707) -
              DENDRO_696 * grad_0_gt2[pp] - DENDRO_697 * grad_1_gt2[pp] -
              DENDRO_698 * grad_2_gt2[pp] + 2.0 * grad_0_Gt0[pp] * gt2[pp] +
              2.0 * grad_0_Gt1[pp] * gt4[pp] + 2.0 * grad_0_Gt2[pp] * gt5[pp] +
              2.0 * grad_2_Gt0[pp] * gt0[pp] + 2.0 * grad_2_Gt1[pp] * gt1[pp] +
              2.0 * grad_2_Gt2[pp] * gt2[pp]) +
         DENDRO_508 * DENDRO_665 - DENDRO_510 * DENDRO_678 -
         DENDRO_512 * DENDRO_677 + DENDRO_516 * DENDRO_679 -
         12 * grad2_0_2_alpha[pp]) +
    DENDRO_700 * grad_0_beta0[pp] + DENDRO_700 * grad_2_beta2[pp] -
    alpha[pp] * (-At2[pp] * K[pp] + DENDRO_31 * DENDRO_702 +
                 DENDRO_37 * DENDRO_704 + DENDRO_40 * DENDRO_705) +
    beta0[pp] * grad_0_At2[pp] + beta1[pp] * grad_1_At2[pp] +
    beta2[pp] * grad_2_At2[pp];
//--
At_rhs11[pp] =
    -At3[pp] * DENDRO_10 + At3[pp] * DENDRO_11 - At3[pp] * DENDRO_7 +
    DENDRO_17 * grad_1_beta0[pp] +
    DENDRO_667 *
        (DENDRO_263 *
             (6.0 * DENDRO_103 * DENDRO_212 * grad_1_gt3[pp] -
              DENDRO_114 * DENDRO_236 * DENDRO_401 + DENDRO_127 * DENDRO_405 +
              DENDRO_135 * (DENDRO_389 - 1.0 * DENDRO_721) -
              DENDRO_135 * (-1.0 * DENDRO_227 * DENDRO_91 + DENDRO_384) -
              DENDRO_135 * (DENDRO_236 * DENDRO_385 + DENDRO_387) +
              DENDRO_141 * (-DENDRO_725 + DENDRO_99 * grad_1_gt3[pp]) +
              DENDRO_141 *
                  (-DENDRO_127 * grad_2_gt3[pp] + DENDRO_244 * DENDRO_86) +
              DENDRO_141 *
                  (DENDRO_227 * grad_1_gt0[pp] - DENDRO_88 * grad_0_gt3[pp]) +
              DENDRO_156 * (DENDRO_393 + DENDRO_684) +
              DENDRO_156 *
                  (-1.0 * DENDRO_725 + 0.25 * DENDRO_99 * grad_1_gt3[pp]) +
              DENDRO_156 * (DENDRO_244 * DENDRO_91 - DENDRO_685) +
              DENDRO_180 * (DENDRO_380 + DENDRO_382) +
              DENDRO_180 * (DENDRO_132 * DENDRO_220 + DENDRO_378 * DENDRO_88) +
              DENDRO_180 * (DENDRO_236 * DENDRO_378 + DENDRO_722) +
              DENDRO_180 * (DENDRO_48 * DENDRO_723 + DENDRO_724) -
              DENDRO_191 * (DENDRO_187 * DENDRO_227 + DENDRO_188 * DENDRO_213 +
                            DENDRO_190 * DENDRO_244 + DENDRO_367) +
              DENDRO_220 * DENDRO_403 + DENDRO_311 * (DENDRO_388 - DENDRO_721) +
              DENDRO_311 *
                  (DENDRO_150 * DENDRO_227 - DENDRO_220 * grad_0_gt3[pp]) +
              DENDRO_311 * (-DENDRO_236 * grad_2_gt3[pp] + DENDRO_399) -
              DENDRO_365 * DENDRO_699 + DENDRO_404 * DENDRO_88 + DENDRO_412 -
              DENDRO_696 * grad_0_gt3[pp] - DENDRO_697 * grad_1_gt3[pp] -
              DENDRO_698 * grad_2_gt3[pp]) -
         DENDRO_366 * DENDRO_56 +
         DENDRO_63 * (DENDRO_244 - DENDRO_58 * DENDRO_719) +
         DENDRO_666 * gt3[pp] +
         DENDRO_720 * (DENDRO_227 - DENDRO_70 * DENDRO_719) -
         12 * grad2_1_1_alpha[pp]) +
    DENDRO_716 * grad_1_beta2[pp] -
    alpha[pp] * (-At3[pp] * K[pp] + DENDRO_31 * DENDRO_674 +
                 DENDRO_672 * DENDRO_717 + DENDRO_676 * DENDRO_718) +
    beta0[pp] * grad_0_At3[pp] + beta1[pp] * grad_1_At3[pp] +
    beta2[pp] * grad_2_At3[pp];
//--
At_rhs12[pp] =
    At1[pp] * grad_2_beta0[pp] + At2[pp] * grad_1_beta0[pp] +
    At3[pp] * grad_2_beta1[pp] - At4[pp] * DENDRO_10 +
    At5[pp] * grad_1_beta2[pp] +
    DENDRO_667 *
        (DENDRO_263 *
             (-DENDRO_104 * (-DENDRO_236 * DENDRO_295 + DENDRO_594) -
              DENDRO_104 *
                  (-DENDRO_100 * DENDRO_220 + DENDRO_550 + DENDRO_595) -
              DENDRO_104 * (0.25 * DENDRO_197 * grad_1_gt3[pp] -
                            DENDRO_210 * DENDRO_212 - DENDRO_211 * DENDRO_212) -
              DENDRO_114 * (DENDRO_203 * DENDRO_240 - 0.5 * DENDRO_727) -
              DENDRO_114 *
                  (-DENDRO_109 * DENDRO_220 + 0.5 * DENDRO_172 * DENDRO_220 -
                   DENDRO_223 * DENDRO_378) +
              DENDRO_135 * (DENDRO_203 * DENDRO_212 + DENDRO_599) -
              DENDRO_135 * (DENDRO_165 * DENDRO_223 + DENDRO_220 * DENDRO_468 +
                            DENDRO_598) +
              DENDRO_135 * (-DENDRO_210 * DENDRO_240 +
                            0.5 * DENDRO_244 * grad_2_gt5[pp] - DENDRO_730) +
              DENDRO_135 * (-DENDRO_220 * DENDRO_483 - DENDRO_220 * DENDRO_575 +
                            0.5 * DENDRO_227 * grad_0_gt5[pp]) +
              DENDRO_135 * (DENDRO_236 * DENDRO_334 + DENDRO_562 - DENDRO_730) +
              DENDRO_156 * (-DENDRO_108 * DENDRO_295 - DENDRO_731) +
              DENDRO_156 * (-DENDRO_212 * DENDRO_318 + DENDRO_601) +
              DENDRO_156 * (-DENDRO_236 * DENDRO_575 - DENDRO_731) +
              DENDRO_156 *
                  (-DENDRO_100 * DENDRO_93 + 0.5 * DENDRO_227 * grad_2_gt0[pp] -
                   0.25 * DENDRO_711) +
              DENDRO_156 * (-DENDRO_178 * DENDRO_220 + DENDRO_227 * DENDRO_51 -
                            DENDRO_724) +
              DENDRO_156 *
                  (-DENDRO_212 * DENDRO_378 + DENDRO_444 + DENDRO_572) +
              DENDRO_160 * (DENDRO_680 + DENDRO_681 + DENDRO_707) -
              DENDRO_180 * (0.5 * DENDRO_108 * DENDRO_203 - DENDRO_729) +
              DENDRO_180 * (DENDRO_240 * DENDRO_378 + DENDRO_729) -
              DENDRO_180 * (-DENDRO_182 * DENDRO_220 - DENDRO_223 * DENDRO_48 -
                            DENDRO_714) -
              DENDRO_191 * (DENDRO_197 * DENDRO_494 + DENDRO_221 * DENDRO_492 +
                            DENDRO_237 * DENDRO_498 + DENDRO_580) +
              DENDRO_298 * (DENDRO_127 * DENDRO_171 + DENDRO_689) +
              DENDRO_298 * (DENDRO_137 + DENDRO_686 + DENDRO_715) +
              DENDRO_311 *
                  (DENDRO_197 * grad_2_gt3[pp] - DENDRO_212 * grad_1_gt5[pp]) -
              DENDRO_586 * DENDRO_699 + DENDRO_587 + DENDRO_607 -
              DENDRO_696 * grad_0_gt4[pp] - DENDRO_697 * grad_1_gt4[pp] -
              DENDRO_698 * grad_2_gt4[pp]) -
         DENDRO_30 * DENDRO_592 * DENDRO_678 + DENDRO_586 * DENDRO_665 +
         DENDRO_590 * DENDRO_679 - DENDRO_591 * DENDRO_677 -
         12 * grad2_1_2_alpha[pp]) +
    DENDRO_726 * grad_1_beta1[pp] + DENDRO_726 * grad_2_beta2[pp] -
    alpha[pp] * (-At4[pp] * K[pp] + DENDRO_31 * DENDRO_704 +
                 DENDRO_702 * DENDRO_717 + DENDRO_705 * DENDRO_718) +
    beta0[pp] * grad_0_At4[pp] + beta1[pp] * grad_1_At4[pp] +
    beta2[pp] * grad_2_At4[pp];
//--
At_rhs22[pp] =
    -At5[pp] * DENDRO_10 + At5[pp] * DENDRO_14 - At5[pp] * DENDRO_9 +
    DENDRO_18 * grad_2_beta0[pp] +
    DENDRO_667 *
        (DENDRO_263 *
             (3.0 * DENDRO_108 * DENDRO_34 * DENDRO_83 * grad_0_gt5[pp] +
              DENDRO_114 * DENDRO_223 * DENDRO_303 -
              DENDRO_135 * (0.25 * DENDRO_727 + 1.0 * DENDRO_733) -
              DENDRO_135 * (-1.0 * DENDRO_224 * DENDRO_86 + DENDRO_337) -
              DENDRO_156 * (DENDRO_108 * DENDRO_200 + DENDRO_728) -
              DENDRO_156 * (DENDRO_136 * DENDRO_220 + DENDRO_318 * DENDRO_93) -
              DENDRO_156 * (DENDRO_171 * DENDRO_236 + DENDRO_708) -
              DENDRO_156 * (DENDRO_51 * DENDRO_723 + 0.25 * DENDRO_682) -
              DENDRO_180 * (-DENDRO_223 * DENDRO_341 - DENDRO_713) -
              DENDRO_191 * (DENDRO_187 * DENDRO_224 + DENDRO_188 * DENDRO_204 +
                            DENDRO_190 * DENDRO_241 + DENDRO_342) +
              4 * DENDRO_20 * DENDRO_319 * DENDRO_83 +
              4 * DENDRO_20 * DENDRO_329 * DENDRO_83 +
              2.0 * DENDRO_22 * DENDRO_310 * DENDRO_83 +
              4 * DENDRO_22 * DENDRO_339 * DENDRO_83 +
              4 * DENDRO_220 * DENDRO_25 * DENDRO_302 * DENDRO_83 +
              3.0 * DENDRO_236 * DENDRO_25 * DENDRO_83 * grad_1_gt5[pp] +
              6.0 * DENDRO_240 * DENDRO_38 * DENDRO_83 * grad_2_gt5[pp] +
              2.0 * DENDRO_25 * DENDRO_30 * grad2_1_1_gt5[pp] -
              DENDRO_272 * DENDRO_699 - DENDRO_278 - DENDRO_279 - DENDRO_280 -
              DENDRO_296 + 4 * DENDRO_297 * DENDRO_34 * DENDRO_83 * DENDRO_93 +
              4 * DENDRO_30 * DENDRO_32 * grad2_0_2_gt5[pp] +
              2.0 * DENDRO_30 * DENDRO_34 * grad2_0_0_gt5[pp] +
              2.0 * DENDRO_30 * DENDRO_38 * grad2_2_2_gt5[pp] - DENDRO_301 -
              DENDRO_304 - DENDRO_311 * (DENDRO_727 + DENDRO_733) -
              DENDRO_311 *
                  (DENDRO_150 * DENDRO_223 + DENDRO_220 * grad_0_gt5[pp]) +
              4 * DENDRO_32 * DENDRO_83 *
                  (0.25 * DENDRO_706 + 1.0 * DENDRO_732) +
              2.0 * DENDRO_32 * DENDRO_83 * (DENDRO_706 + DENDRO_732) +
              4 * DENDRO_32 * DENDRO_83 *
                  (-1.0 * DENDRO_204 * DENDRO_86 + DENDRO_335) +
              2.0 * DENDRO_32 * DENDRO_83 *
                  (DENDRO_223 * grad_2_gt0[pp] + DENDRO_93 * grad_0_gt5[pp]) -
              DENDRO_325 - DENDRO_696 * grad_0_gt5[pp] -
              DENDRO_697 * grad_1_gt5[pp] - DENDRO_698 * grad_2_gt5[pp] +
              4 * grad_2_Gt0[pp] * gt2[pp] + 4 * grad_2_Gt1[pp] * gt4[pp] +
              4 * grad_2_Gt2[pp] * gt5[pp]) +
         DENDRO_265 * DENDRO_57 - DENDRO_268 * DENDRO_720 -
         DENDRO_274 * DENDRO_62 + DENDRO_666 * gt5[pp] -
         12 * grad2_2_2_alpha[pp]) +
    DENDRO_716 * grad_2_beta1[pp] -
    alpha[pp] * (At5[pp] * DENDRO_36 * DENDRO_705 - At5[pp] * K[pp] +
                 DENDRO_40 * DENDRO_704 + DENDRO_702 * DENDRO_718) +
    beta0[pp] * grad_0_At5[pp] + beta1[pp] * grad_1_At5[pp] +
    beta2[pp] * grad_2_At5[pp];
//--
K_rhs[pp] =
    DENDRO_252 * DENDRO_30 * chi[pp] *
        (0.5 * DENDRO_741 * (DENDRO_636 + DENDRO_655 * DENDRO_739) +
         DENDRO_743 *
             (DENDRO_30 * DENDRO_634 + DENDRO_30 * DENDRO_96 +
              DENDRO_30 * DENDRO_97 -
              DENDRO_44 * (-DENDRO_651 * DENDRO_737 + grad_0_chi[pp])) +
         DENDRO_745 *
             (DENDRO_30 * DENDRO_632 + DENDRO_30 * DENDRO_633 +
              DENDRO_30 * DENDRO_84 -
              DENDRO_44 * (-DENDRO_651 * DENDRO_734 + grad_1_chi[pp])) -
         grad2_0_1_alpha[pp]) +
    DENDRO_255 * DENDRO_30 * chi[pp] *
        (0.5 * DENDRO_736 * (DENDRO_44 * DENDRO_734 * gt4[pp] + DENDRO_581) +
         DENDRO_743 * (DENDRO_30 * DENDRO_582 -
                       DENDRO_44 * (-DENDRO_586 * DENDRO_737 + grad_2_chi[pp]) +
                       DENDRO_589) +
         DENDRO_744 *
             (DENDRO_30 * DENDRO_583 + DENDRO_30 * DENDRO_584 +
              DENDRO_30 * DENDRO_585 -
              DENDRO_44 * (-DENDRO_586 * DENDRO_739 + grad_1_chi[pp])) -
         grad2_1_2_alpha[pp]) -
    DENDRO_283 * chi[pp] *
        (DENDRO_738 * (DENDRO_430 + DENDRO_45 * DENDRO_742) +
         DENDRO_741 * (DENDRO_431 + DENDRO_45 * DENDRO_740) -
         grad2_0_0_alpha[pp] +
         grad_0_alpha[pp] *
             (DENDRO_30 * DENDRO_428 + DENDRO_30 * DENDRO_429 -
              DENDRO_44 * (DENDRO_68 - 0.5 * DENDRO_69 * DENDRO_734) +
              DENDRO_67)) -
    DENDRO_285 * chi[pp] *
        (DENDRO_736 * (DENDRO_369 + DENDRO_719 * DENDRO_734) +
         DENDRO_741 * (DENDRO_373 + DENDRO_719 * DENDRO_739) -
         grad2_1_1_alpha[pp] +
         grad_1_alpha[pp] *
             (DENDRO_30 * DENDRO_370 + DENDRO_30 * DENDRO_371 +
              DENDRO_30 * DENDRO_372 -
              DENDRO_44 * (DENDRO_364 - DENDRO_365 * DENDRO_742))) -
    DENDRO_287 * chi[pp] *
        (DENDRO_736 * (DENDRO_347 + DENDRO_734 * DENDRO_735) +
         DENDRO_738 * (DENDRO_349 + DENDRO_735 * DENDRO_737) -
         grad2_2_2_alpha[pp] +
         grad_2_alpha[pp] *
             (DENDRO_30 * DENDRO_351 + DENDRO_30 * DENDRO_352 +
              DENDRO_30 * DENDRO_353 -
              DENDRO_44 * (DENDRO_271 - DENDRO_272 * DENDRO_740))) -
    DENDRO_746 * chi[pp] *
        (0.5 * DENDRO_738 * (DENDRO_493 + DENDRO_515 * DENDRO_737) +
         DENDRO_744 *
             (DENDRO_30 * DENDRO_495 + DENDRO_30 * DENDRO_496 +
              DENDRO_30 * DENDRO_497 -
              DENDRO_44 * (-DENDRO_508 * DENDRO_739 + grad_0_chi[pp])) +
         DENDRO_745 *
             (DENDRO_30 * DENDRO_490 + DENDRO_30 * DENDRO_491 +
              DENDRO_30 * DENDRO_92 -
              DENDRO_44 * (-DENDRO_508 * DENDRO_734 + grad_2_chi[pp])) -
         grad2_0_2_alpha[pp]) +
    DENDRO_765 *
        (At0[pp] * DENDRO_750 * DENDRO_751 + At1[pp] * DENDRO_760 * DENDRO_762 +
         At2[pp] * DENDRO_759 * DENDRO_760 + At3[pp] * DENDRO_751 * DENDRO_754 +
         At4[pp] * DENDRO_760 * DENDRO_764 + At5[pp] * DENDRO_751 * DENDRO_757 +
         pow(K[pp], 2)) +
    beta0[pp] * grad_0_K[pp] + beta1[pp] * grad_1_K[pp] +
    beta2[pp] * grad_2_K[pp];
//--
Gt_rhs0[pp] = -DENDRO_780;
//--
Gt_rhs1[pp] =
    DENDRO_116 * DENDRO_759 * DENDRO_799 - 2.0 / 3.0 * DENDRO_16 * DENDRO_216 -
    DENDRO_197 * DENDRO_763 * DENDRO_799 +
    DENDRO_204 * DENDRO_757 * DENDRO_798 -
    DENDRO_212 * DENDRO_754 * DENDRO_798 + DENDRO_216 * grad_1_beta1[pp] +
    DENDRO_231 * grad_0_beta1[pp] + DENDRO_248 * grad_2_beta1[pp] +
    DENDRO_54 * DENDRO_750 * DENDRO_798 + DENDRO_761 * DENDRO_770 -
    DENDRO_761 * DENDRO_799 * DENDRO_99 + DENDRO_763 * DENDRO_774 +
    DENDRO_773 * DENDRO_793 -
    DENDRO_773 * (-DENDRO_761 * DENDRO_787 + DENDRO_786) -
    DENDRO_773 * (-DENDRO_763 * DENDRO_789 + DENDRO_788) - DENDRO_781 -
    DENDRO_782 - DENDRO_783 - DENDRO_784 - DENDRO_785 - DENDRO_791 -
    DENDRO_792 + (7.0 / 3.0) * DENDRO_794 * grad2_0_1_beta1[pp] +
    DENDRO_795 * grad2_0_0_beta0[pp] + DENDRO_795 * grad2_0_2_beta2[pp] +
    DENDRO_796 * DENDRO_797 + (7.0 / 3.0) * DENDRO_796 * grad2_1_2_beta1[pp] +
    (1.0 / 3.0) * DENDRO_796 * grad2_2_2_beta2[pp] + DENDRO_800;
//--
Gt_rhs2[pp] = -DENDRO_801;
//--
B_rhs0[pp] =
    -B0[pp] * eta - DENDRO_780 +
    lambda[2] * (beta0[pp] * grad_0_B0[pp] + beta1[pp] * grad_1_B0[pp] +
                 beta2[pp] * grad_2_B0[pp]) -
    lambda[3] * (beta0[pp] * grad_0_Gt0[pp] + beta1[pp] * grad_1_Gt0[pp] +
                 beta2[pp] * grad_2_Gt0[pp]);
//--
B_rhs1[pp] =
    -B1[pp] * eta + 2.0 * DENDRO_116 * DENDRO_759 * DENDRO_771 * alpha[pp] +
    (2.0 / 3.0) * DENDRO_16 * DENDRO_258 * DENDRO_83 +
    2.0 * DENDRO_197 * DENDRO_764 * DENDRO_771 * alpha[pp] +
    (1.0 / 3.0) * DENDRO_20 * DENDRO_30 * grad2_0_0_beta0[pp] +
    (7.0 / 3.0) * DENDRO_20 * DENDRO_30 * grad2_0_1_beta1[pp] +
    (1.0 / 3.0) * DENDRO_20 * DENDRO_30 * grad2_0_2_beta2[pp] +
    2 * DENDRO_204 * DENDRO_757 * DENDRO_771 * alpha[pp] +
    2 * DENDRO_213 * DENDRO_754 * DENDRO_771 * alpha[pp] +
    (1.0 / 3.0) * DENDRO_22 * DENDRO_30 * grad2_0_2_beta0[pp] +
    (7.0 / 3.0) * DENDRO_22 * DENDRO_30 * grad2_1_2_beta1[pp] +
    (1.0 / 3.0) * DENDRO_22 * DENDRO_30 * grad2_2_2_beta2[pp] -
    DENDRO_354 * grad_1_beta1[pp] - DENDRO_356 * grad_0_beta1[pp] -
    DENDRO_358 * grad_2_beta1[pp] +
    2 * DENDRO_54 * DENDRO_750 * DENDRO_771 * alpha[pp] -
    DENDRO_762 * DENDRO_770 +
    2.0 * DENDRO_762 * DENDRO_771 * DENDRO_99 * alpha[pp] -
    DENDRO_764 * DENDRO_774 + DENDRO_773 * DENDRO_793 -
    DENDRO_773 * (DENDRO_762 * DENDRO_787 + DENDRO_786) -
    DENDRO_773 * (DENDRO_764 * DENDRO_789 + DENDRO_788) - DENDRO_781 -
    DENDRO_782 - DENDRO_783 - DENDRO_784 - DENDRO_785 - DENDRO_791 -
    DENDRO_792 - DENDRO_800 * lambda[3] + beta0[pp] * grad_0_Gt1[pp] +
    beta1[pp] * grad_1_Gt1[pp] + beta2[pp] * grad_2_Gt1[pp] +
    lambda[2] * (beta0[pp] * grad_0_B1[pp] + beta1[pp] * grad_1_B1[pp] +
                 beta2[pp] * grad_2_B1[pp]);
//--
B_rhs2[pp] =
    -B2[pp] * eta - DENDRO_801 +
    lambda[2] * (beta0[pp] * grad_0_B2[pp] + beta1[pp] * grad_1_B2[pp] +
                 beta2[pp] * grad_2_B2[pp]) -
    lambda[3] * (beta0[pp] * grad_0_Gt2[pp] + beta1[pp] * grad_1_Gt2[pp] +
                 beta2[pp] * grad_2_Gt2[pp]);
// Dendro: reduced ops: 4302
// Dendro: }}}
