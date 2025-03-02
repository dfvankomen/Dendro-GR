// Codgen: generating unstage version
// Codgen: using standard gauge
// Codgen: using eta const damping
//  Dendro: {{{
//  Dendro: original ops: 623308
//  Dendro: printing temp variables
const double DENDRO_0 = 2 * alpha[pp];
const double DENDRO_1 =
    (3.0 / 4.0) * alpha[pp] * lambda_f[1] + (3.0 / 4.0) * lambda_f[0];
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
const double DENDRO_24 = gt0[pp] * gt5[pp];
const double DENDRO_25 = pow(gt2[pp], 2);
const double DENDRO_26 = DENDRO_24 - DENDRO_25;
const double DENDRO_27 = At1[pp] * DENDRO_26;
const double DENDRO_28 = pow(gt4[pp], 2);
const double DENDRO_29 = pow(gt1[pp], 2);
const double DENDRO_30 = gt3[pp] * gt5[pp];
const double DENDRO_31 = -DENDRO_19 * DENDRO_3 + DENDRO_25 * gt3[pp] +
                         DENDRO_28 * gt0[pp] + DENDRO_29 * gt5[pp] -
                         DENDRO_30 * gt0[pp];
const double DENDRO_32 = 1.0 / DENDRO_31;
const double DENDRO_33 = DENDRO_17 * DENDRO_32;
const double DENDRO_34 = At1[pp] * DENDRO_20;
const double DENDRO_35 = gt1[pp] * gt4[pp];
const double DENDRO_36 = gt2[pp] * gt3[pp];
const double DENDRO_37 = DENDRO_35 - DENDRO_36;
const double DENDRO_38 = At2[pp] * DENDRO_37;
const double DENDRO_39 = -DENDRO_38;
const double DENDRO_40 = -DENDRO_28 + DENDRO_30;
const double DENDRO_41 = At0[pp] * DENDRO_40;
const double DENDRO_42 = 2 * DENDRO_32;
const double DENDRO_43 = At0[pp] * DENDRO_42;
const double DENDRO_44 = At1[pp] * DENDRO_22;
const double DENDRO_45 = gt0[pp] * gt3[pp];
const double DENDRO_46 = -DENDRO_29 + DENDRO_45;
const double DENDRO_47 = At2[pp] * DENDRO_46;
const double DENDRO_48 = DENDRO_18 * DENDRO_32;
const double DENDRO_49 =
    DENDRO_20 * grad_0_chi[pp] + DENDRO_22 * grad_2_chi[pp];
const double DENDRO_50  = -DENDRO_26 * grad_1_chi[pp] + DENDRO_49;
const double DENDRO_51  = 0.5 * DENDRO_50;
const double DENDRO_52  = 1.0 / chi[pp];
const double DENDRO_53  = DENDRO_52 * gt0[pp];
const double DENDRO_54  = 1.0 * grad_0_gt1[pp];
const double DENDRO_55  = 0.5 * grad_1_gt0[pp];
const double DENDRO_56  = DENDRO_54 - DENDRO_55;
const double DENDRO_57  = 1.0 * grad_0_gt2[pp];
const double DENDRO_58  = 0.5 * grad_2_gt0[pp];
const double DENDRO_59  = DENDRO_57 - DENDRO_58;
const double DENDRO_60  = 0.5 * grad_0_gt0[pp];
const double DENDRO_61  = DENDRO_20 * DENDRO_60 + DENDRO_22 * DENDRO_59;
const double DENDRO_62  = -DENDRO_26 * DENDRO_56 + DENDRO_61;
const double DENDRO_63  = DENDRO_51 * DENDRO_53 + DENDRO_62;
const double DENDRO_64  = 12 * grad_1_alpha[pp];
const double DENDRO_65  = DENDRO_32 * DENDRO_64;
const double DENDRO_66  = DENDRO_22 * grad_1_chi[pp];
const double DENDRO_67  = DENDRO_37 * grad_0_chi[pp];
const double DENDRO_68  = DENDRO_46 * grad_2_chi[pp];
const double DENDRO_69  = DENDRO_66 - DENDRO_67 - DENDRO_68;
const double DENDRO_70  = 0.5 * DENDRO_69;
const double DENDRO_71  = DENDRO_53 * DENDRO_70;
const double DENDRO_72  = DENDRO_37 * DENDRO_60;
const double DENDRO_73  = DENDRO_22 * DENDRO_56;
const double DENDRO_74  = DENDRO_46 * DENDRO_59;
const double DENDRO_75  = DENDRO_72 - DENDRO_73 + DENDRO_74;
const double DENDRO_76  = 12 * grad_2_alpha[pp];
const double DENDRO_77  = DENDRO_32 * DENDRO_76;
const double DENDRO_78  = DENDRO_40 * DENDRO_60;
const double DENDRO_79  = DENDRO_32 * DENDRO_78;
const double DENDRO_80  = DENDRO_37 * DENDRO_59;
const double DENDRO_81  = DENDRO_32 * DENDRO_80;
const double DENDRO_82  = DENDRO_20 * DENDRO_56;
const double DENDRO_83  = DENDRO_32 * DENDRO_82;
const double DENDRO_84  = 1.0 * grad_0_chi[pp];
const double DENDRO_85  = DENDRO_32 * gt0[pp];
const double DENDRO_86  = DENDRO_20 * grad_1_chi[pp];
const double DENDRO_87  = DENDRO_37 * grad_2_chi[pp];
const double DENDRO_88  = DENDRO_40 * grad_0_chi[pp];
const double DENDRO_89  = DENDRO_86 - DENDRO_87 - DENDRO_88;
const double DENDRO_90  = 0.5 * DENDRO_89;
const double DENDRO_91  = DENDRO_52 * (DENDRO_84 - DENDRO_85 * DENDRO_90);
const double DENDRO_92  = 12 * grad_0_alpha[pp];
const double DENDRO_93  = -grad2_0_0_chi[pp];
const double DENDRO_94  = -DENDRO_78 - DENDRO_80 + DENDRO_82;
const double DENDRO_95  = DENDRO_32 * grad_0_chi[pp];
const double DENDRO_96  = DENDRO_32 * grad_1_chi[pp];
const double DENDRO_97  = -DENDRO_72 + DENDRO_73 - DENDRO_74;
const double DENDRO_98  = DENDRO_32 * grad_2_chi[pp];
const double DENDRO_99  = 2 * DENDRO_52;
const double DENDRO_100 = 0.5 * grad_0_gt3[pp];
const double DENDRO_101 = 1.0 * grad_1_gt1[pp];
const double DENDRO_102 = DENDRO_100 - DENDRO_101;
const double DENDRO_103 = grad_0_gt4[pp] - grad_1_gt2[pp] + grad_2_gt1[pp];
const double DENDRO_104 =
    DENDRO_20 * grad_2_gt0[pp] + DENDRO_22 * grad_0_gt5[pp];
const double DENDRO_105 = -DENDRO_103 * DENDRO_26 + DENDRO_104;
const double DENDRO_106 = DENDRO_102 * DENDRO_105;
const double DENDRO_107 = DENDRO_26 * grad_0_gt3[pp];
const double DENDRO_108 = DENDRO_20 * grad_1_gt0[pp];
const double DENDRO_109 = grad_0_gt4[pp] + grad_1_gt2[pp] - grad_2_gt1[pp];
const double DENDRO_110 = DENDRO_109 * DENDRO_22;
const double DENDRO_111 = DENDRO_108 + DENDRO_110;
const double DENDRO_112 = -DENDRO_107 + DENDRO_111;
const double DENDRO_113 = DENDRO_103 * DENDRO_112;
const double DENDRO_114 = 0.25 * DENDRO_113;
const double DENDRO_115 = pow(DENDRO_31, -2);
const double DENDRO_116 = DENDRO_115 * DENDRO_22;
const double DENDRO_117 = 4 * DENDRO_116;
const double DENDRO_118 = 0.5 * grad_0_gt5[pp];
const double DENDRO_119 = 1.0 * grad_2_gt2[pp];
const double DENDRO_120 = -DENDRO_119;
const double DENDRO_121 = DENDRO_118 + DENDRO_120;
const double DENDRO_122 = DENDRO_22 * grad_0_gt3[pp];
const double DENDRO_123 = DENDRO_37 * grad_1_gt0[pp];
const double DENDRO_124 = DENDRO_109 * DENDRO_46;
const double DENDRO_125 = DENDRO_122 - DENDRO_123 - DENDRO_124;
const double DENDRO_126 = DENDRO_121 * DENDRO_125;
const double DENDRO_127 = DENDRO_37 * grad_2_gt0[pp];
const double DENDRO_128 = DENDRO_46 * grad_0_gt5[pp];
const double DENDRO_129 = DENDRO_103 * DENDRO_22;
const double DENDRO_130 = -DENDRO_127 - DENDRO_128 + DENDRO_129;
const double DENDRO_131 = DENDRO_109 * DENDRO_130;
const double DENDRO_132 = 0.25 * DENDRO_131;
const double DENDRO_133 = -DENDRO_132;
const double DENDRO_134 = -DENDRO_122 + DENDRO_123 + DENDRO_124;
const double DENDRO_135 = 0.25 * grad_0_gt5[pp];
const double DENDRO_136 = DENDRO_134 * DENDRO_135;
const double DENDRO_137 = -grad_0_gt4[pp] + grad_1_gt2[pp] + grad_2_gt1[pp];
const double DENDRO_138 = DENDRO_127 + DENDRO_128 - DENDRO_129;
const double DENDRO_139 = 0.5 * DENDRO_138;
const double DENDRO_140 = DENDRO_37 * grad_0_gt5[pp];
const double DENDRO_141 = DENDRO_40 * grad_2_gt0[pp];
const double DENDRO_142 = DENDRO_103 * DENDRO_20;
const double DENDRO_143 = DENDRO_140 + DENDRO_141 - DENDRO_142;
const double DENDRO_144 = 0.25 * grad_1_gt0[pp];
const double DENDRO_145 = DENDRO_143 * DENDRO_144;
const double DENDRO_146 = DENDRO_20 * grad_0_gt3[pp];
const double DENDRO_147 = DENDRO_40 * grad_1_gt0[pp];
const double DENDRO_148 = DENDRO_109 * DENDRO_37;
const double DENDRO_149 = -DENDRO_146 + DENDRO_147 + DENDRO_148;
const double DENDRO_150 = 0.25 * grad_2_gt0[pp];
const double DENDRO_151 = DENDRO_149 * DENDRO_150;
const double DENDRO_152 = 2 * DENDRO_121;
const double DENDRO_153 = DENDRO_115 * DENDRO_37;
const double DENDRO_154 = 4 * DENDRO_153;
const double DENDRO_155 = DENDRO_149 * grad_0_gt0[pp];
const double DENDRO_156 = 0.25 * DENDRO_155;
const double DENDRO_157 = DENDRO_78 + DENDRO_80 - DENDRO_82;
const double DENDRO_158 = DENDRO_157 * grad_1_gt0[pp];
const double DENDRO_159 = DENDRO_115 * DENDRO_20;
const double DENDRO_160 = 4 * DENDRO_159;
const double DENDRO_161 = 0.5 * DENDRO_59;
const double DENDRO_162 = DENDRO_143 * grad_0_gt0[pp];
const double DENDRO_163 = 0.25 * DENDRO_162;
const double DENDRO_164 = DENDRO_157 * grad_2_gt0[pp];
const double DENDRO_165 = 2.0 * DENDRO_153;
const double DENDRO_166 = DENDRO_138 * grad_2_gt0[pp];
const double DENDRO_167 = DENDRO_75 * grad_0_gt5[pp];
const double DENDRO_168 = -DENDRO_140 - DENDRO_141 + DENDRO_142;
const double DENDRO_169 = DENDRO_168 * DENDRO_37;
const double DENDRO_170 = DENDRO_20 * grad_2_gt3[pp];
const double DENDRO_171 = DENDRO_37 * grad_1_gt5[pp];
const double DENDRO_172 = DENDRO_137 * DENDRO_40;
const double DENDRO_173 = DENDRO_170 - DENDRO_171 - DENDRO_172;
const double DENDRO_174 = DENDRO_173 * DENDRO_22;
const double DENDRO_175 = DENDRO_146 - DENDRO_147 - DENDRO_148;
const double DENDRO_176 = DENDRO_175 * DENDRO_20;
const double DENDRO_177 = 0.5 * grad_2_gt5[pp];
const double DENDRO_178 = DENDRO_177 * DENDRO_37;
const double DENDRO_179 = DENDRO_121 * DENDRO_40;
const double DENDRO_180 = 0.5 * grad_1_gt5[pp];
const double DENDRO_181 = 1.0 * grad_2_gt4[pp];
const double DENDRO_182 = -DENDRO_181;
const double DENDRO_183 = DENDRO_180 + DENDRO_182;
const double DENDRO_184 = DENDRO_183 * DENDRO_20;
const double DENDRO_185 = -DENDRO_178 + DENDRO_179 - DENDRO_184;
const double DENDRO_186 = DENDRO_185 * DENDRO_46;
const double DENDRO_187 = 0.5 * grad_1_gt3[pp];
const double DENDRO_188 = DENDRO_187 * DENDRO_20;
const double DENDRO_189 = 1.0 * grad_1_gt4[pp];
const double DENDRO_190 = 0.5 * grad_2_gt3[pp];
const double DENDRO_191 = DENDRO_189 - DENDRO_190;
const double DENDRO_192 =
    DENDRO_102 * DENDRO_40 + DENDRO_188 - DENDRO_191 * DENDRO_37;
const double DENDRO_193 = DENDRO_192 * DENDRO_26;
const double DENDRO_194 = DENDRO_40 * DENDRO_94;
const double DENDRO_195 =
    DENDRO_115 * (DENDRO_169 - 1.0 * DENDRO_174 - 1.0 * DENDRO_176 +
                  DENDRO_186 + DENDRO_193 + DENDRO_194);
const double DENDRO_196 = 2.0 * grad_0_gt0[pp];
const double DENDRO_197 = DENDRO_105 * DENDRO_37;
const double DENDRO_198 = DENDRO_26 * grad_2_gt3[pp];
const double DENDRO_199 = DENDRO_22 * grad_1_gt5[pp];
const double DENDRO_200 = DENDRO_137 * DENDRO_20;
const double DENDRO_201 = DENDRO_199 + DENDRO_200;
const double DENDRO_202 = -DENDRO_198 + DENDRO_201;
const double DENDRO_203 = DENDRO_202 * DENDRO_22;
const double DENDRO_204 = DENDRO_112 * DENDRO_20;
const double DENDRO_205 = DENDRO_177 * DENDRO_22;
const double DENDRO_206 =
    -DENDRO_121 * DENDRO_20 + DENDRO_183 * DENDRO_26 + DENDRO_205;
const double DENDRO_207 = DENDRO_206 * DENDRO_46;
const double DENDRO_208 = DENDRO_187 * DENDRO_26;
const double DENDRO_209 = DENDRO_191 * DENDRO_22;
const double DENDRO_210 = DENDRO_102 * DENDRO_20;
const double DENDRO_211 = -DENDRO_208 + DENDRO_209 - DENDRO_210;
const double DENDRO_212 = DENDRO_211 * DENDRO_26;
const double DENDRO_213 = DENDRO_40 * DENDRO_62;
const double DENDRO_214 =
    DENDRO_115 * (DENDRO_197 - 1.0 * DENDRO_203 - 1.0 * DENDRO_204 +
                  DENDRO_207 + DENDRO_212 + DENDRO_213);
const double DENDRO_215 = 2.0 * grad_1_gt0[pp];
const double DENDRO_216 = DENDRO_130 * DENDRO_37;
const double DENDRO_217 = DENDRO_22 * grad_2_gt3[pp];
const double DENDRO_218 = DENDRO_46 * grad_1_gt5[pp];
const double DENDRO_219 = DENDRO_137 * DENDRO_37;
const double DENDRO_220 = DENDRO_217 - DENDRO_218 - DENDRO_219;
const double DENDRO_221 = DENDRO_22 * DENDRO_220;
const double DENDRO_222 = DENDRO_125 * DENDRO_20;
const double DENDRO_223 = DENDRO_177 * DENDRO_46;
const double DENDRO_224 = DENDRO_121 * DENDRO_37;
const double DENDRO_225 = DENDRO_183 * DENDRO_22;
const double DENDRO_226 = -DENDRO_223 + DENDRO_224 - DENDRO_225;
const double DENDRO_227 = DENDRO_226 * DENDRO_46;
const double DENDRO_228 = DENDRO_187 * DENDRO_22;
const double DENDRO_229 =
    DENDRO_102 * DENDRO_37 - DENDRO_191 * DENDRO_46 + DENDRO_228;
const double DENDRO_230 = DENDRO_229 * DENDRO_26;
const double DENDRO_231 = DENDRO_40 * DENDRO_97;
const double DENDRO_232 =
    DENDRO_115 * (DENDRO_216 - 1.0 * DENDRO_221 - 1.0 * DENDRO_222 +
                  DENDRO_227 + DENDRO_230 + DENDRO_231);
const double DENDRO_233 = 2.0 * grad_2_gt0[pp];
const double DENDRO_234 = 2.0 * DENDRO_159;
const double DENDRO_235 = DENDRO_134 * grad_2_gt0[pp];
const double DENDRO_236 = DENDRO_109 * DENDRO_75;
const double DENDRO_237 = 3 * DENDRO_52;
const double DENDRO_238 = grad_0_chi[pp] * grad_1_chi[pp];
const double DENDRO_239 = 2 * DENDRO_20;
const double DENDRO_240 =
    DENDRO_239 * (-DENDRO_237 * DENDRO_238 + 2 * grad2_0_1_chi[pp]);
const double DENDRO_241 = DENDRO_237 * grad_2_chi[pp];
const double DENDRO_242 = 2 * DENDRO_37;
const double DENDRO_243 =
    DENDRO_242 * (-DENDRO_241 * grad_0_chi[pp] + 2 * grad2_0_2_chi[pp]);
const double DENDRO_244 = 2 * DENDRO_22;
const double DENDRO_245 =
    DENDRO_244 * (-DENDRO_241 * grad_1_chi[pp] + 2 * grad2_1_2_chi[pp]);
const double DENDRO_246 = pow(grad_0_chi[pp], 2);
const double DENDRO_247 =
    DENDRO_40 * (-DENDRO_237 * DENDRO_246 + 2 * grad2_0_0_chi[pp]);
const double DENDRO_248 = pow(grad_1_chi[pp], 2);
const double DENDRO_249 =
    DENDRO_26 * (-DENDRO_237 * DENDRO_248 + 2 * grad2_1_1_chi[pp]);
const double DENDRO_250 = pow(grad_2_chi[pp], 2);
const double DENDRO_251 =
    DENDRO_46 * (-DENDRO_237 * DENDRO_250 + 2 * grad2_2_2_chi[pp]);
const double DENDRO_252 = -1.0 * DENDRO_197 + DENDRO_203 + DENDRO_204 -
                          DENDRO_207 - DENDRO_212 - DENDRO_213;
const double DENDRO_253 = 2 * DENDRO_252 * DENDRO_96;
const double DENDRO_254 = -1.0 * DENDRO_169 + DENDRO_174 + DENDRO_176 -
                          DENDRO_186 - DENDRO_193 - DENDRO_194;
const double DENDRO_255 = 2 * DENDRO_254 * DENDRO_95;
const double DENDRO_256 = -1.0 * DENDRO_216 + DENDRO_221 + DENDRO_222 -
                          DENDRO_227 - DENDRO_230 - DENDRO_231;
const double DENDRO_257 = 2 * DENDRO_256 * DENDRO_98;
const double DENDRO_258 = -DENDRO_240 + DENDRO_243 - DENDRO_245 + DENDRO_247 +
                          DENDRO_249 + DENDRO_251 + DENDRO_253 + DENDRO_255 +
                          DENDRO_257;
const double DENDRO_259 = DENDRO_52 * DENDRO_85;
const double DENDRO_260 = DENDRO_115 * DENDRO_46;
const double DENDRO_261 = 4 * DENDRO_260;
const double DENDRO_262 = 0.25 * grad_0_gt4[pp];
const double DENDRO_263 = -DENDRO_262;
const double DENDRO_264 = 0.75 * grad_1_gt2[pp];
const double DENDRO_265 = 0.25 * grad_2_gt1[pp];
const double DENDRO_266 = DENDRO_115 * DENDRO_26;
const double DENDRO_267 = 4 * DENDRO_266;
const double DENDRO_268 = DENDRO_267 * (DENDRO_263 + DENDRO_264 + DENDRO_265);
const double DENDRO_269 = DENDRO_57 + DENDRO_58;
const double DENDRO_270 = DENDRO_115 * DENDRO_40;
const double DENDRO_271 = 4 * DENDRO_270;
const double DENDRO_272 = DENDRO_149 * grad_1_gt0[pp];
const double DENDRO_273 = 3.0 * DENDRO_266;
const double DENDRO_274 = DENDRO_143 * grad_2_gt0[pp];
const double DENDRO_275 = 3.0 * DENDRO_260;
const double DENDRO_276 = 6.0 * grad_0_gt0[pp];
const double DENDRO_277 = pow(chi[pp], -2);
const double DENDRO_278 = 4 * gt1[pp];
const double DENDRO_279 = 4 * gt2[pp];
const double DENDRO_280 = 0.5 * DENDRO_56;
const double DENDRO_281 = DENDRO_105 * DENDRO_280;
const double DENDRO_282 = DENDRO_32 * DENDRO_37;
const double DENDRO_283 = 4 * DENDRO_282;
const double DENDRO_284 = 0.25 * grad_0_gt3[pp];
const double DENDRO_285 = DENDRO_105 * DENDRO_284;
const double DENDRO_286 = 0.5 * DENDRO_137;
const double DENDRO_287 = DENDRO_112 * DENDRO_280;
const double DENDRO_288 = DENDRO_32 * DENDRO_40;
const double DENDRO_289 = 2.0 * DENDRO_288;
const double DENDRO_290 = DENDRO_26 * DENDRO_32;
const double DENDRO_291 = 2.0 * DENDRO_290;
const double DENDRO_292 = DENDRO_32 * DENDRO_46;
const double DENDRO_293 = 2.0 * DENDRO_292;
const double DENDRO_294 = DENDRO_112 * grad_1_gt0[pp];
const double DENDRO_295 = DENDRO_62 * grad_0_gt3[pp];
const double DENDRO_296 = DENDRO_105 * grad_1_gt0[pp];
const double DENDRO_297 = DENDRO_103 * DENDRO_62;
const double DENDRO_298 = 4.0 * DENDRO_32;
const double DENDRO_299 = DENDRO_20 * DENDRO_298;
const double DENDRO_300 = DENDRO_22 * DENDRO_298;
const double DENDRO_301 = 0.25 * grad_1_gt2[pp];
const double DENDRO_302 = 0.75 * grad_2_gt1[pp];
const double DENDRO_303 = 4 * DENDRO_115;
const double DENDRO_304 =
    -DENDRO_105 * DENDRO_261 * (DENDRO_263 + DENDRO_301 + DENDRO_302) -
    DENDRO_112 * DENDRO_267 * (DENDRO_101 - DENDRO_284) +
    DENDRO_117 * (DENDRO_112 * DENDRO_286 + DENDRO_285) -
    DENDRO_154 * (DENDRO_137 * DENDRO_62 + DENDRO_281) +
    DENDRO_160 * (-2 * DENDRO_102 * DENDRO_62 + DENDRO_287) -
    DENDRO_165 * (DENDRO_296 + DENDRO_297) -
    DENDRO_213 * DENDRO_303 * (DENDRO_54 + DENDRO_55) +
    DENDRO_234 * (DENDRO_294 + DENDRO_295) - DENDRO_246 * DENDRO_277 +
    DENDRO_278 * grad_0_Gt1[pp] + DENDRO_279 * grad_0_Gt2[pp] +
    DENDRO_283 * grad2_0_2_gt0[pp] + DENDRO_289 * grad2_0_0_gt0[pp] +
    DENDRO_291 * grad2_1_1_gt0[pp] + DENDRO_293 * grad2_2_2_gt0[pp] -
    DENDRO_299 * grad2_0_1_gt0[pp] - DENDRO_300 * grad2_1_2_gt0[pp] +
    4 * grad_0_Gt0[pp] * gt0[pp];
const double DENDRO_305 = 3 * alpha[pp];
const double DENDRO_306 = DENDRO_52 * gt5[pp];
const double DENDRO_307 = DENDRO_206 + DENDRO_306 * DENDRO_51;
const double DENDRO_308 = 4 * grad_1_alpha[pp];
const double DENDRO_309 = DENDRO_308 * DENDRO_32;
const double DENDRO_310 = DENDRO_306 * DENDRO_90;
const double DENDRO_311 = 4 * grad_0_alpha[pp];
const double DENDRO_312 = DENDRO_311 * DENDRO_32;
const double DENDRO_313 = DENDRO_223 * DENDRO_32;
const double DENDRO_314 = DENDRO_224 * DENDRO_32;
const double DENDRO_315 = DENDRO_225 * DENDRO_32;
const double DENDRO_316 = 1.0 * grad_2_chi[pp];
const double DENDRO_317 = DENDRO_32 * gt5[pp];
const double DENDRO_318 = DENDRO_52 * (DENDRO_316 - DENDRO_317 * DENDRO_70);
const double DENDRO_319 = 4 * grad_2_alpha[pp];
const double DENDRO_320 = -grad2_2_2_chi[pp];
const double DENDRO_321 = -DENDRO_35 + DENDRO_36;
const double DENDRO_322 = -DENDRO_118;
const double DENDRO_323 = DENDRO_119 + DENDRO_322;
const double DENDRO_324 = DENDRO_28 - DENDRO_30;
const double DENDRO_325 = -DENDRO_180 + DENDRO_181;
const double DENDRO_326 =
    DENDRO_177 * DENDRO_321 + DENDRO_20 * DENDRO_325 + DENDRO_323 * DENDRO_324;
const double DENDRO_327 = -DENDRO_24 + DENDRO_25;
const double DENDRO_328 =
    DENDRO_20 * DENDRO_323 + DENDRO_205 + DENDRO_325 * DENDRO_327;
const double DENDRO_329 = DENDRO_29 - DENDRO_45;
const double DENDRO_330 = DENDRO_177 * DENDRO_329;
const double DENDRO_331 = DENDRO_321 * DENDRO_323;
const double DENDRO_332 = DENDRO_22 * DENDRO_325;
const double DENDRO_333 = 0.5 * DENDRO_183;
const double DENDRO_334 = DENDRO_105 * DENDRO_333;
const double DENDRO_335 = -DENDRO_334;
const double DENDRO_336 = DENDRO_109 * DENDRO_206;
const double DENDRO_337 = 0.5 * DENDRO_121;
const double DENDRO_338 = -DENDRO_168 * DENDRO_337;
const double DENDRO_339 = 2 * DENDRO_59;
const double DENDRO_340 = DENDRO_130 * grad_2_gt5[pp];
const double DENDRO_341 = 0.25 * DENDRO_340;
const double DENDRO_342 = DENDRO_226 * grad_0_gt5[pp];
const double DENDRO_343 = DENDRO_173 * DENDRO_337;
const double DENDRO_344 = -DENDRO_343;
const double DENDRO_345 = DENDRO_109 * DENDRO_185;
const double DENDRO_346 = DENDRO_220 * grad_2_gt5[pp];
const double DENDRO_347 = 0.25 * DENDRO_346;
const double DENDRO_348 = DENDRO_226 * grad_1_gt5[pp];
const double DENDRO_349 = DENDRO_173 * DENDRO_59;
const double DENDRO_350 = DENDRO_137 * DENDRO_168;
const double DENDRO_351 = 0.25 * DENDRO_350;
const double DENDRO_352 = DENDRO_135 * DENDRO_220;
const double DENDRO_353 = 0.25 * grad_1_gt5[pp];
const double DENDRO_354 = DENDRO_130 * DENDRO_353;
const double DENDRO_355 = DENDRO_150 * DENDRO_173;
const double DENDRO_356 = 0.5 * DENDRO_109;
const double DENDRO_357 = DENDRO_115 * DENDRO_254;
const double DENDRO_358 = 2.0 * DENDRO_357;
const double DENDRO_359 = DENDRO_115 * DENDRO_252;
const double DENDRO_360 = 2.0 * DENDRO_359;
const double DENDRO_361 = DENDRO_115 * DENDRO_256;
const double DENDRO_362 = 2.0 * DENDRO_361;
const double DENDRO_363 = DENDRO_137 * DENDRO_185;
const double DENDRO_364 = 2.0 * DENDRO_116;
const double DENDRO_365 = DENDRO_185 * grad_2_gt0[pp];
const double DENDRO_366 = DENDRO_240 - DENDRO_243 + DENDRO_245 - DENDRO_247 -
                          DENDRO_249 - DENDRO_251 - DENDRO_253 - DENDRO_255 -
                          DENDRO_257;
const double DENDRO_367 = DENDRO_366 * DENDRO_52;
const double DENDRO_368 = -DENDRO_265;
const double DENDRO_369 = DENDRO_267 * (DENDRO_262 + DENDRO_264 + DENDRO_368);
const double DENDRO_370 = DENDRO_271 * (-DENDRO_150 + DENDRO_57);
const double DENDRO_371 = DENDRO_130 * grad_0_gt5[pp];
const double DENDRO_372 = 3.0 * DENDRO_270;
const double DENDRO_373 = DENDRO_220 * grad_1_gt5[pp];
const double DENDRO_374 = 6.0 * DENDRO_115;
const double DENDRO_375 = 4 * gt4[pp];
const double DENDRO_376 = -DENDRO_202 * DENDRO_333;
const double DENDRO_377 = DENDRO_105 * DENDRO_191;
const double DENDRO_378 = DENDRO_103 * DENDRO_202;
const double DENDRO_379 = 0.25 * DENDRO_378;
const double DENDRO_380 = 0.25 * grad_2_gt3[pp];
const double DENDRO_381 = DENDRO_105 * DENDRO_380;
const double DENDRO_382 = DENDRO_202 * grad_1_gt5[pp];
const double DENDRO_383 = DENDRO_206 * grad_2_gt3[pp];
const double DENDRO_384 = DENDRO_105 * grad_1_gt5[pp];
const double DENDRO_385 = DENDRO_103 * DENDRO_206;
const double DENDRO_386 = 0.75 * grad_0_gt4[pp];
const double DENDRO_387 =
    -DENDRO_105 * DENDRO_271 * (DENDRO_301 + DENDRO_368 + DENDRO_386) +
    DENDRO_117 * (2 * DENDRO_191 * DENDRO_206 + DENDRO_376) +
    DENDRO_160 * (DENDRO_377 + DENDRO_379) +
    DENDRO_160 * (DENDRO_202 * DENDRO_356 + DENDRO_381) -
    DENDRO_165 * (DENDRO_384 + DENDRO_385) -
    DENDRO_202 * DENDRO_267 * (DENDRO_189 - DENDRO_380) -
    DENDRO_207 * DENDRO_303 * (DENDRO_180 + DENDRO_181) -
    DENDRO_250 * DENDRO_277 + DENDRO_279 * grad_2_Gt0[pp] +
    DENDRO_283 * grad2_0_2_gt5[pp] + DENDRO_289 * grad2_0_0_gt5[pp] +
    DENDRO_291 * grad2_1_1_gt5[pp] + DENDRO_293 * grad2_2_2_gt5[pp] -
    DENDRO_299 * grad2_0_1_gt5[pp] - DENDRO_300 * grad2_1_2_gt5[pp] +
    DENDRO_364 * (DENDRO_382 + DENDRO_383) + DENDRO_375 * grad_2_Gt1[pp] +
    4 * grad_2_Gt2[pp] * gt5[pp];
const double DENDRO_388 = DENDRO_52 * gt3[pp];
const double DENDRO_389 = DENDRO_319 * DENDRO_32;
const double DENDRO_390 = DENDRO_208 * DENDRO_32;
const double DENDRO_391 = DENDRO_209 * DENDRO_32;
const double DENDRO_392 = DENDRO_210 * DENDRO_32;
const double DENDRO_393 = 1.0 * grad_1_chi[pp];
const double DENDRO_394 = DENDRO_32 * gt3[pp];
const double DENDRO_395 = DENDRO_52 * (DENDRO_393 - DENDRO_394 * DENDRO_51);
const double DENDRO_396 = -grad2_1_1_chi[pp];
const double DENDRO_397 = -DENDRO_100 + DENDRO_101;
const double DENDRO_398 =
    DENDRO_188 + DENDRO_191 * DENDRO_321 + DENDRO_324 * DENDRO_397;
const double DENDRO_399 = DENDRO_187 * DENDRO_327;
const double DENDRO_400 = DENDRO_20 * DENDRO_397;
const double DENDRO_401 =
    DENDRO_191 * DENDRO_329 + DENDRO_228 + DENDRO_321 * DENDRO_397;
const double DENDRO_402 = DENDRO_173 * DENDRO_56;
const double DENDRO_403 = DENDRO_137 * DENDRO_175;
const double DENDRO_404 = 0.25 * DENDRO_403;
const double DENDRO_405 = DENDRO_144 * DENDRO_173;
const double DENDRO_406 = 0.5 * DENDRO_103;
const double DENDRO_407 = DENDRO_125 * DENDRO_353;
const double DENDRO_408 = DENDRO_109 * DENDRO_220;
const double DENDRO_409 = 0.25 * DENDRO_408;
const double DENDRO_410 = DENDRO_125 * DENDRO_183;
const double DENDRO_411 = 0.5 * DENDRO_102;
const double DENDRO_412 = DENDRO_173 * DENDRO_411;
const double DENDRO_413 = -DENDRO_412;
const double DENDRO_414 = DENDRO_103 * DENDRO_192;
const double DENDRO_415 = 0.5 * DENDRO_191;
const double DENDRO_416 = DENDRO_220 * DENDRO_415;
const double DENDRO_417 = 2 * DENDRO_183 * DENDRO_229;
const double DENDRO_418 = DENDRO_202 * grad_1_gt3[pp];
const double DENDRO_419 = 0.25 * DENDRO_418;
const double DENDRO_420 = DENDRO_211 * grad_2_gt3[pp];
const double DENDRO_421 = DENDRO_125 * DENDRO_415;
const double DENDRO_422 = DENDRO_103 * DENDRO_229;
const double DENDRO_423 = -DENDRO_175 * DENDRO_411;
const double DENDRO_424 = 2 * DENDRO_192 * DENDRO_56;
const double DENDRO_425 = DENDRO_112 * grad_1_gt3[pp];
const double DENDRO_426 = 0.25 * DENDRO_425;
const double DENDRO_427 = DENDRO_211 * grad_0_gt3[pp];
const double DENDRO_428 = DENDRO_192 * grad_1_gt0[pp];
const double DENDRO_429 = DENDRO_137 * DENDRO_192;
const double DENDRO_430 = DENDRO_229 * grad_1_gt5[pp];
const double DENDRO_431 = DENDRO_109 * DENDRO_229;
const double DENDRO_432 = -DENDRO_301;
const double DENDRO_433 = DENDRO_261 * (DENDRO_262 + DENDRO_302 + DENDRO_432);
const double DENDRO_434 = DENDRO_271 * (-DENDRO_144 + DENDRO_54);
const double DENDRO_435 = DENDRO_271 * (DENDRO_265 + DENDRO_386 + DENDRO_432);
const double DENDRO_436 = DENDRO_202 * DENDRO_284;
const double DENDRO_437 = DENDRO_112 * DENDRO_380;
const double DENDRO_438 = DENDRO_112 * grad_0_gt3[pp];
const double DENDRO_439 = DENDRO_202 * grad_2_gt3[pp];
const double DENDRO_440 =
    -DENDRO_154 * (DENDRO_100 * DENDRO_202 + DENDRO_437) -
    DENDRO_154 * (DENDRO_112 * DENDRO_190 + DENDRO_436) -
    DENDRO_193 * DENDRO_303 * (DENDRO_100 + DENDRO_101) -
    DENDRO_230 * DENDRO_303 * (DENDRO_189 + DENDRO_190) -
    DENDRO_248 * DENDRO_277 - DENDRO_275 * DENDRO_439 +
    DENDRO_278 * grad_1_Gt0[pp] + DENDRO_283 * grad2_0_2_gt3[pp] +
    DENDRO_289 * grad2_0_0_gt3[pp] + DENDRO_291 * grad2_1_1_gt3[pp] +
    DENDRO_293 * grad2_2_2_gt3[pp] - DENDRO_299 * grad2_0_1_gt3[pp] -
    DENDRO_300 * grad2_1_2_gt3[pp] - DENDRO_372 * DENDRO_438 +
    DENDRO_375 * grad_1_Gt2[pp] + 4 * grad_1_Gt1[pp] * gt3[pp];
const double DENDRO_441 = DENDRO_324 * DENDRO_60;
const double DENDRO_442 = DENDRO_321 * DENDRO_59;
const double DENDRO_443 = DENDRO_327 * DENDRO_56 + DENDRO_61;
const double DENDRO_444 =
    DENDRO_321 * DENDRO_60 + DENDRO_329 * DENDRO_59 + DENDRO_73;
const double DENDRO_445 = DENDRO_130 * DENDRO_161;
const double DENDRO_446 = DENDRO_168 * grad_0_gt0[pp];
const double DENDRO_447 = 0.25 * DENDRO_446;
const double DENDRO_448 = DENDRO_94 * grad_2_gt0[pp];
const double DENDRO_449 = DENDRO_125 * DENDRO_135;
const double DENDRO_450 = DENDRO_144 * DENDRO_168;
const double DENDRO_451 = DENDRO_150 * DENDRO_175;
const double DENDRO_452 = DENDRO_125 * DENDRO_161;
const double DENDRO_453 = DENDRO_175 * grad_0_gt0[pp];
const double DENDRO_454 = 0.25 * DENDRO_453;
const double DENDRO_455 = DENDRO_94 * grad_1_gt0[pp];
const double DENDRO_456 = DENDRO_109 * DENDRO_97;
const double DENDRO_457 = DENDRO_97 * grad_0_gt5[pp];
const double DENDRO_458 = DENDRO_175 * grad_1_gt0[pp];
const double DENDRO_459 = DENDRO_168 * grad_2_gt0[pp];
const double DENDRO_460 = DENDRO_202 * grad_1_gt0[pp];
const double DENDRO_461 = DENDRO_105 * grad_0_gt3[pp];
const double DENDRO_462 = DENDRO_109 * DENDRO_112;
const double DENDRO_463 = DENDRO_461 + DENDRO_462;
const double DENDRO_464 = DENDRO_220 * grad_2_gt0[pp];
const double DENDRO_465 = DENDRO_125 * grad_0_gt5[pp];
const double DENDRO_466 = DENDRO_131 + DENDRO_465;
const double DENDRO_467 = DENDRO_135 * DENDRO_168;
const double DENDRO_468 = -DENDRO_121 * DENDRO_226;
const double DENDRO_469 = DENDRO_206 * DENDRO_286 + 0.25 * DENDRO_384;
const double DENDRO_470 = DENDRO_109 * DENDRO_175;
const double DENDRO_471 = 0.25 * DENDRO_470;
const double DENDRO_472 = DENDRO_112 * DENDRO_415;
const double DENDRO_473 = -DENDRO_202 * DENDRO_411;
const double DENDRO_474 = DENDRO_472 + DENDRO_473;
const double DENDRO_475 = DENDRO_59 * DENDRO_94;
const double DENDRO_476 = 0.25 * DENDRO_296 + DENDRO_356 * DENDRO_62;
const double DENDRO_477 = DENDRO_130 * DENDRO_150;
const double DENDRO_478 = 0.25 * DENDRO_105;
const double DENDRO_479 = DENDRO_137 * DENDRO_478 + DENDRO_180 * DENDRO_62;
const double DENDRO_480 = -DENDRO_130 * DENDRO_337;
const double DENDRO_481 = DENDRO_177 * DENDRO_97;
const double DENDRO_482 = DENDRO_150 * DENDRO_168;
const double DENDRO_483 = DENDRO_185 * DENDRO_60;
const double DENDRO_484 = DENDRO_161 * DENDRO_168 + DENDRO_483;
const double DENDRO_485 = DENDRO_103 * DENDRO_478;
const double DENDRO_486 = DENDRO_109 * DENDRO_478 + DENDRO_206 * DENDRO_55;
const double DENDRO_487 = DENDRO_102 * DENDRO_206;
const double DENDRO_488 = 0.5 * DENDRO_377;
const double DENDRO_489 = -DENDRO_487 + DENDRO_488;
const double DENDRO_490 = DENDRO_185 * DENDRO_55;
const double DENDRO_491 = DENDRO_135 * DENDRO_175 + DENDRO_490;
const double DENDRO_492 = DENDRO_226 * DENDRO_286;
const double DENDRO_493 = DENDRO_352 + DENDRO_354;
const double DENDRO_494 = DENDRO_137 * DENDRO_202;
const double DENDRO_495 = 0.25 * DENDRO_494;
const double DENDRO_496 = DENDRO_100 * DENDRO_206;
const double DENDRO_497 = DENDRO_112 * DENDRO_353 + DENDRO_496;
const double DENDRO_498 = 0.25 * DENDRO_109;
const double DENDRO_499 = DENDRO_168 * DENDRO_498 + DENDRO_490;
const double DENDRO_500 = DENDRO_220 * DENDRO_337;
const double DENDRO_501 = -DENDRO_500;
const double DENDRO_502 = 0.25 * DENDRO_125;
const double DENDRO_503 = DENDRO_502 * grad_2_gt5[pp];
const double DENDRO_504 = DENDRO_226 * DENDRO_356 + DENDRO_503;
const double DENDRO_505 = DENDRO_202 * DENDRO_280;
const double DENDRO_506 = -0.5 * DENDRO_106 + DENDRO_191 * DENDRO_62;
const double DENDRO_507 = 0.25 * DENDRO_173;
const double DENDRO_508 = DENDRO_507 * grad_0_gt0[pp];
const double DENDRO_509 = DENDRO_161 * DENDRO_175;
const double DENDRO_510 = DENDRO_508 + DENDRO_509;
const double DENDRO_511 = DENDRO_356 * DENDRO_94 + DENDRO_450;
const double DENDRO_512 = 0.25 * DENDRO_137;
const double DENDRO_513 = DENDRO_130 * DENDRO_512;
const double DENDRO_514 = DENDRO_161 * DENDRO_220;
const double DENDRO_515 = DENDRO_180 * DENDRO_97;
const double DENDRO_516 = DENDRO_137 * DENDRO_220;
const double DENDRO_517 = DENDRO_125 * grad_1_gt5[pp];
const double DENDRO_518 = DENDRO_408 + DENDRO_517;
const double DENDRO_519 = 1.0 * DENDRO_266;
const double DENDRO_520 = -grad2_0_2_chi[pp];
const double DENDRO_521 = DENDRO_321 * grad_0_gt5[pp];
const double DENDRO_522 = DENDRO_324 * grad_2_gt0[pp];
const double DENDRO_523 = 0.5 * DENDRO_95;
const double DENDRO_524 = DENDRO_103 * DENDRO_327 + DENDRO_104;
const double DENDRO_525 = 0.5 * DENDRO_96;
const double DENDRO_526 = DENDRO_329 * grad_0_gt5[pp];
const double DENDRO_527 = DENDRO_321 * grad_2_gt0[pp];
const double DENDRO_528 = 0.5 * DENDRO_98;
const double DENDRO_529 = 2.0 * gt2[pp];
const double DENDRO_530 = DENDRO_529 * grad_0_Gt0[pp];
const double DENDRO_531 = 2.0 * gt4[pp];
const double DENDRO_532 = DENDRO_531 * grad_0_Gt1[pp];
const double DENDRO_533 = 2.0 * gt5[pp];
const double DENDRO_534 = DENDRO_533 * grad_0_Gt2[pp];
const double DENDRO_535 = 2.0 * grad_2_Gt0[pp];
const double DENDRO_536 = DENDRO_535 * gt0[pp];
const double DENDRO_537 = 2.0 * grad_2_Gt1[pp];
const double DENDRO_538 = DENDRO_537 * gt1[pp];
const double DENDRO_539 = DENDRO_529 * grad_2_Gt2[pp];
const double DENDRO_540 = DENDRO_277 * grad_2_chi[pp];
const double DENDRO_541 = -DENDRO_540 * grad_0_chi[pp];
const double DENDRO_542 = DENDRO_283 * grad2_0_2_gt2[pp];
const double DENDRO_543 = DENDRO_289 * grad2_0_0_gt2[pp];
const double DENDRO_544 = DENDRO_291 * grad2_1_1_gt2[pp];
const double DENDRO_545 = DENDRO_293 * grad2_2_2_gt2[pp];
const double DENDRO_546 = -DENDRO_299 * grad2_0_1_gt2[pp];
const double DENDRO_547 = -DENDRO_300 * grad2_1_2_gt2[pp];
const double DENDRO_548 = DENDRO_32 * gt2[pp];
const double DENDRO_549 =
    DENDRO_358 * grad_0_gt2[pp] + DENDRO_360 * grad_1_gt2[pp] +
    DENDRO_362 * grad_2_gt2[pp] + DENDRO_367 * DENDRO_548 + DENDRO_530 +
    DENDRO_532 + DENDRO_534 + DENDRO_536 + DENDRO_538 + DENDRO_539 +
    DENDRO_541 + DENDRO_542 + DENDRO_543 + DENDRO_544 + DENDRO_545 +
    DENDRO_546 + DENDRO_547 -
    DENDRO_99 *
        (DENDRO_520 + DENDRO_523 * (DENDRO_142 + DENDRO_521 + DENDRO_522) +
         DENDRO_524 * DENDRO_525 +
         DENDRO_528 * (DENDRO_129 + DENDRO_526 + DENDRO_527));
const double DENDRO_550 = DENDRO_140 * DENDRO_32;
const double DENDRO_551 = DENDRO_141 * DENDRO_32;
const double DENDRO_552 = DENDRO_142 * DENDRO_32;
const double DENDRO_553 =
    DENDRO_52 * (-DENDRO_548 * DENDRO_89 + grad_2_chi[pp]);
const double DENDRO_554 = 2.0 * grad_0_alpha[pp];
const double DENDRO_555 = DENDRO_127 * DENDRO_32;
const double DENDRO_556 = DENDRO_128 * DENDRO_32;
const double DENDRO_557 = DENDRO_129 * DENDRO_32;
const double DENDRO_558 =
    DENDRO_52 * (-DENDRO_548 * DENDRO_69 + grad_0_chi[pp]);
const double DENDRO_559 = 2.0 * grad_2_alpha[pp];
const double DENDRO_560 = 2.0 * grad_1_alpha[pp];
const double DENDRO_561 = DENDRO_52 * gt2[pp];
const double DENDRO_562 = DENDRO_32 * (DENDRO_105 + DENDRO_50 * DENDRO_561);
const double DENDRO_563 =
    DENDRO_554 * (-DENDRO_550 - DENDRO_551 + DENDRO_552 - DENDRO_553) +
    DENDRO_559 * (-DENDRO_555 - DENDRO_556 + DENDRO_557 - DENDRO_558) +
    DENDRO_560 * DENDRO_562 - 4 * grad2_0_2_alpha[pp];
const double DENDRO_564 = DENDRO_173 * grad_2_gt0[pp];
const double DENDRO_565 = DENDRO_175 * grad_0_gt5[pp] + DENDRO_564;
const double DENDRO_566 = 0.25 * DENDRO_516;
const double DENDRO_567 = DENDRO_130 * DENDRO_135 + DENDRO_481;
const double DENDRO_568 = DENDRO_112 * DENDRO_333;
const double DENDRO_569 = -DENDRO_568;
const double DENDRO_570 = DENDRO_354 + DENDRO_503;
const double DENDRO_571 = -DENDRO_175 * DENDRO_337;
const double DENDRO_572 = DENDRO_286 * DENDRO_94;
const double DENDRO_573 = DENDRO_150 * DENDRO_220;
const double DENDRO_574 = DENDRO_449 + DENDRO_515;
const double DENDRO_575 = DENDRO_173 * grad_1_gt0[pp];
const double DENDRO_576 = DENDRO_470 + DENDRO_575;
const double DENDRO_577 = DENDRO_105 * grad_2_gt3[pp];
const double DENDRO_578 = DENDRO_112 * grad_1_gt5[pp] + DENDRO_577;
const double DENDRO_579 = 0.25 * DENDRO_462;
const double DENDRO_580 = DENDRO_190 * DENDRO_62;
const double DENDRO_581 = DENDRO_144 * DENDRO_202 + DENDRO_580;
const double DENDRO_582 = DENDRO_116 * (DENDRO_494 + DENDRO_578) -
                          DENDRO_154 * (DENDRO_479 + DENDRO_485) -
                          DENDRO_154 * (-DENDRO_183 * DENDRO_62 + DENDRO_486) +
                          DENDRO_160 * (DENDRO_114 + DENDRO_506) +
                          DENDRO_160 * (DENDRO_579 + DENDRO_581) -
                          DENDRO_261 * (DENDRO_335 + DENDRO_469) -
                          DENDRO_267 * (DENDRO_437 + DENDRO_474) -
                          DENDRO_271 * (0.5 * DENDRO_297 + DENDRO_476);
const double DENDRO_583 = DENDRO_168 * grad_0_gt3[pp];
const double DENDRO_584 = DENDRO_130 * grad_2_gt3[pp];
const double DENDRO_585 = 0.25 * DENDRO_382;
const double DENDRO_586 = -DENDRO_183 * DENDRO_226;
const double DENDRO_587 = DENDRO_135 * DENDRO_173 + DENDRO_185 * DENDRO_406;
const double DENDRO_588 = DENDRO_191 * DENDRO_211;
const double DENDRO_589 = DENDRO_192 * DENDRO_356;
const double DENDRO_590 = DENDRO_173 * DENDRO_284 + DENDRO_589;
const double DENDRO_591 = DENDRO_220 * DENDRO_380;
const double DENDRO_592 = DENDRO_168 * DENDRO_280;
const double DENDRO_593 = DENDRO_509 + DENDRO_592;
const double DENDRO_594 = DENDRO_185 * DENDRO_56 + 0.5 * DENDRO_349;
const double DENDRO_595 = DENDRO_103 * DENDRO_168;
const double DENDRO_596 = 0.25 * DENDRO_595;
const double DENDRO_597 = DENDRO_226 * DENDRO_406;
const double DENDRO_598 = DENDRO_202 * DENDRO_498 + DENDRO_496;
const double DENDRO_599 = DENDRO_130 * DENDRO_333;
const double DENDRO_600 = -DENDRO_599;
const double DENDRO_601 = DENDRO_177 * DENDRO_229;
const double DENDRO_602 = -DENDRO_220 * DENDRO_333 + DENDRO_601;
const double DENDRO_603 = DENDRO_118 * DENDRO_192;
const double DENDRO_604 = DENDRO_103 * DENDRO_507 + DENDRO_603;
const double DENDRO_605 = DENDRO_187 * DENDRO_206;
const double DENDRO_606 = DENDRO_202 * DENDRO_380 + DENDRO_605;
const double DENDRO_607 = DENDRO_202 * DENDRO_415;
const double DENDRO_608 = DENDRO_137 * DENDRO_507;
const double DENDRO_609 = DENDRO_100 * DENDRO_185 + DENDRO_173 * DENDRO_498;
const double DENDRO_610 = -DENDRO_168 * DENDRO_411;
const double DENDRO_611 = DENDRO_192 * DENDRO_59;
const double DENDRO_612 = 0.5 * DENDRO_402 + DENDRO_611;
const double DENDRO_613 = DENDRO_478 * grad_1_gt3[pp];
const double DENDRO_614 = DENDRO_436 + DENDRO_613;
const double DENDRO_615 = DENDRO_211 * DENDRO_356;
const double DENDRO_616 = DENDRO_130 * DENDRO_415;
const double DENDRO_617 = 0.25 * DENDRO_103;
const double DENDRO_618 = DENDRO_118 * DENDRO_229;
const double DENDRO_619 = DENDRO_220 * DENDRO_617 + DENDRO_618;
const double DENDRO_620 = DENDRO_103 * DENDRO_130;
const double DENDRO_621 = 1.0 * DENDRO_270;
const double DENDRO_622 = -grad2_1_2_chi[pp];
const double DENDRO_623 =
    DENDRO_137 * DENDRO_324 + DENDRO_170 + DENDRO_321 * grad_1_gt5[pp];
const double DENDRO_624 = DENDRO_327 * grad_2_gt3[pp];
const double DENDRO_625 = DENDRO_329 * grad_1_gt5[pp];
const double DENDRO_626 = DENDRO_137 * DENDRO_321;
const double DENDRO_627 = DENDRO_32 * gt4[pp];
const double DENDRO_628 =
    DENDRO_283 * grad2_0_2_gt4[pp] + DENDRO_289 * grad2_0_0_gt4[pp] +
    DENDRO_291 * grad2_1_1_gt4[pp] + DENDRO_293 * grad2_2_2_gt4[pp] -
    DENDRO_299 * grad2_0_1_gt4[pp] - DENDRO_300 * grad2_1_2_gt4[pp] +
    DENDRO_529 * grad_1_Gt0[pp] + DENDRO_531 * grad_1_Gt1[pp] +
    DENDRO_531 * grad_2_Gt2[pp] + DENDRO_533 * grad_1_Gt2[pp] +
    DENDRO_535 * gt1[pp] + DENDRO_537 * gt3[pp] - DENDRO_540 * grad_1_chi[pp];
const double DENDRO_629 =
    DENDRO_358 * grad_0_gt4[pp] + DENDRO_360 * grad_1_gt4[pp] +
    DENDRO_362 * grad_2_gt4[pp] + DENDRO_367 * DENDRO_627 + DENDRO_628 -
    DENDRO_99 *
        (DENDRO_523 * DENDRO_623 + DENDRO_525 * (DENDRO_201 + DENDRO_624) +
         DENDRO_528 * (DENDRO_217 + DENDRO_625 + DENDRO_626) + DENDRO_622);
const double DENDRO_630 = DENDRO_199 * DENDRO_32 + DENDRO_200 * DENDRO_32;
const double DENDRO_631 =
    -DENDRO_198 * DENDRO_32 -
    DENDRO_52 * (-DENDRO_50 * DENDRO_627 + grad_2_chi[pp]) + DENDRO_630;
const double DENDRO_632 = DENDRO_217 * DENDRO_32;
const double DENDRO_633 = DENDRO_218 * DENDRO_32;
const double DENDRO_634 = DENDRO_219 * DENDRO_32;
const double DENDRO_635 =
    DENDRO_52 * (-DENDRO_627 * DENDRO_69 + grad_1_chi[pp]);
const double DENDRO_636 = DENDRO_52 * gt4[pp];
const double DENDRO_637 = DENDRO_636 * DENDRO_89;
const double DENDRO_638 =
    DENDRO_32 * DENDRO_554 * (DENDRO_173 + DENDRO_637) +
    DENDRO_559 * (DENDRO_632 - DENDRO_633 - DENDRO_634 - DENDRO_635) +
    DENDRO_560 * DENDRO_631 - 4 * grad2_1_2_alpha[pp];
const double DENDRO_639 = 1.0 * DENDRO_430;
const double DENDRO_640 = 0.5 * DENDRO_429;
const double DENDRO_641 = 0.25 * DENDRO_620;
const double DENDRO_642 = DENDRO_352 + DENDRO_503;
const double DENDRO_643 = DENDRO_121 * DENDRO_192;
const double DENDRO_644 = DENDRO_605 + DENDRO_607;
const double DENDRO_645 = DENDRO_220 * DENDRO_353;
const double DENDRO_646 = DENDRO_436 + DENDRO_437;
const double DENDRO_647 = DENDRO_192 * DENDRO_58;
const double DENDRO_648 = DENDRO_168 * DENDRO_284 + DENDRO_647;
const double DENDRO_649 = DENDRO_211 * DENDRO_406 + DENDRO_613;
const double DENDRO_650 = DENDRO_130 * DENDRO_380 + DENDRO_618;
const double DENDRO_651 = 1.0 * DENDRO_153;
const double DENDRO_652 =
    -DENDRO_154 * (DENDRO_569 + DENDRO_598) -
    DENDRO_261 * (DENDRO_190 * DENDRO_206 + DENDRO_376 + DENDRO_585) -
    DENDRO_621 * (DENDRO_113 + DENDRO_463) -
    DENDRO_651 * (DENDRO_378 + DENDRO_578);
const double DENDRO_653 = DENDRO_501 + DENDRO_600;
const double DENDRO_654 = DENDRO_175 * DENDRO_284;
const double DENDRO_655 = -DENDRO_102 * DENDRO_211;
const double DENDRO_656 = DENDRO_229 * DENDRO_286;
const double DENDRO_657 = DENDRO_125 * DENDRO_380 + DENDRO_656;
const double DENDRO_658 = DENDRO_56 * DENDRO_94;
const double DENDRO_659 = 0.25 * DENDRO_294;
const double DENDRO_660 = DENDRO_125 * DENDRO_150 + DENDRO_406 * DENDRO_97;
const double DENDRO_661 = DENDRO_183 * DENDRO_97;
const double DENDRO_662 = 0.5 * DENDRO_126;
const double DENDRO_663 = -DENDRO_661 - DENDRO_662;
const double DENDRO_664 = DENDRO_508 + DENDRO_592;
const double DENDRO_665 = DENDRO_406 * DENDRO_94 + DENDRO_451;
const double DENDRO_666 = DENDRO_112 * DENDRO_512 + DENDRO_580;
const double DENDRO_667 = DENDRO_121 * DENDRO_229;
const double DENDRO_668 = 0.5 * DENDRO_410;
const double DENDRO_669 = -DENDRO_667 - DENDRO_668;
const double DENDRO_670 = DENDRO_211 * DENDRO_286;
const double DENDRO_671 = DENDRO_175 * DENDRO_617 + DENDRO_647;
const double DENDRO_672 = -DENDRO_112 * DENDRO_411;
const double DENDRO_673 = DENDRO_187 * DENDRO_62;
const double DENDRO_674 = DENDRO_137 * DENDRO_502 + DENDRO_190 * DENDRO_97;
const double DENDRO_675 = DENDRO_144 * DENDRO_175;
const double DENDRO_676 = DENDRO_192 * DENDRO_60;
const double DENDRO_677 = DENDRO_175 * DENDRO_280 + DENDRO_676;
const double DENDRO_678 = DENDRO_109 * DENDRO_502;
const double DENDRO_679 = DENDRO_229 * DENDRO_58;
const double DENDRO_680 = DENDRO_103 * DENDRO_502 + DENDRO_679;
const double DENDRO_681 = 1.0 * DENDRO_260;
const double DENDRO_682 = -grad2_0_1_chi[pp];
const double DENDRO_683 = DENDRO_324 * grad_1_gt0[pp];
const double DENDRO_684 = DENDRO_109 * DENDRO_321;
const double DENDRO_685 = DENDRO_327 * grad_0_gt3[pp];
const double DENDRO_686 =
    DENDRO_109 * DENDRO_329 + DENDRO_122 + DENDRO_321 * grad_1_gt0[pp];
const double DENDRO_687 = 2.0 * gt1[pp];
const double DENDRO_688 = DENDRO_687 * grad_0_Gt0[pp];
const double DENDRO_689 = 2.0 * grad_0_Gt1[pp] * gt3[pp];
const double DENDRO_690 = DENDRO_531 * grad_0_Gt2[pp];
const double DENDRO_691 = 2.0 * grad_1_Gt0[pp] * gt0[pp];
const double DENDRO_692 = DENDRO_687 * grad_1_Gt1[pp];
const double DENDRO_693 = DENDRO_529 * grad_1_Gt2[pp];
const double DENDRO_694 = -DENDRO_238 * DENDRO_277;
const double DENDRO_695 = DENDRO_283 * grad2_0_2_gt1[pp];
const double DENDRO_696 = DENDRO_289 * grad2_0_0_gt1[pp];
const double DENDRO_697 = DENDRO_291 * grad2_1_1_gt1[pp];
const double DENDRO_698 = DENDRO_293 * grad2_2_2_gt1[pp];
const double DENDRO_699 = -DENDRO_299 * grad2_0_1_gt1[pp];
const double DENDRO_700 = -DENDRO_300 * grad2_1_2_gt1[pp];
const double DENDRO_701 = DENDRO_32 * gt1[pp];
const double DENDRO_702 =
    DENDRO_358 * grad_0_gt1[pp] + DENDRO_360 * grad_1_gt1[pp] +
    DENDRO_362 * grad_2_gt1[pp] + DENDRO_367 * DENDRO_701 + DENDRO_688 +
    DENDRO_689 + DENDRO_690 + DENDRO_691 + DENDRO_692 + DENDRO_693 +
    DENDRO_694 + DENDRO_695 + DENDRO_696 + DENDRO_697 + DENDRO_698 +
    DENDRO_699 + DENDRO_700 -
    DENDRO_99 * (DENDRO_523 * (DENDRO_146 + DENDRO_683 + DENDRO_684) +
                 DENDRO_525 * (DENDRO_111 + DENDRO_685) +
                 DENDRO_528 * DENDRO_686 + DENDRO_682);
const double DENDRO_703 = DENDRO_146 * DENDRO_32;
const double DENDRO_704 = DENDRO_147 * DENDRO_32;
const double DENDRO_705 = DENDRO_148 * DENDRO_32;
const double DENDRO_706 =
    DENDRO_52 * (-DENDRO_701 * DENDRO_89 + grad_1_chi[pp]);
const double DENDRO_707 =
    DENDRO_52 * (-DENDRO_50 * DENDRO_701 + grad_0_chi[pp]);
const double DENDRO_708 = DENDRO_107 * DENDRO_32;
const double DENDRO_709 = DENDRO_108 * DENDRO_32;
const double DENDRO_710 = DENDRO_110 * DENDRO_32;
const double DENDRO_711 = DENDRO_709 + DENDRO_710;
const double DENDRO_712 = DENDRO_52 * gt1[pp];
const double DENDRO_713 =
    DENDRO_32 * DENDRO_559 * (DENDRO_125 + DENDRO_69 * DENDRO_712) +
    DENDRO_554 * (DENDRO_703 - DENDRO_704 - DENDRO_705 - DENDRO_706) +
    DENDRO_560 * (-DENDRO_707 - DENDRO_708 + DENDRO_711) -
    4 * grad2_0_1_alpha[pp];
const double DENDRO_714 = 0.5 * DENDRO_425;
const double DENDRO_715 = DENDRO_192 * DENDRO_55;
const double DENDRO_716 = DENDRO_437 + DENDRO_613;
const double DENDRO_717 = DENDRO_112 * DENDRO_284 + DENDRO_673;
const double DENDRO_718 = DENDRO_117 * (DENDRO_473 + DENDRO_716) -
                          DENDRO_154 * (DENDRO_285 + DENDRO_581) -
                          DENDRO_154 * (DENDRO_285 + DENDRO_666) +
                          DENDRO_160 * (DENDRO_672 + DENDRO_717) -
                          DENDRO_261 * (DENDRO_105 * DENDRO_190 + DENDRO_495) -
                          DENDRO_271 * (1.0 * DENDRO_295 + DENDRO_659);
const double DENDRO_719 =
    -DENDRO_20 *
        (DENDRO_713 +
         alpha[pp] *
             (DENDRO_116 * (DENDRO_403 + DENDRO_575 + DENDRO_583) +
              DENDRO_116 * (DENDRO_516 + DENDRO_517 + DENDRO_584) +
              DENDRO_117 * (DENDRO_610 + DENDRO_671) +
              DENDRO_117 * (DENDRO_616 + DENDRO_669) +
              DENDRO_117 * (DENDRO_670 + DENDRO_716) -
              DENDRO_154 * (DENDRO_132 + DENDRO_663) -
              DENDRO_154 * (DENDRO_450 + DENDRO_665) -
              DENDRO_154 * (DENDRO_572 + DENDRO_664) -
              DENDRO_154 * (DENDRO_515 + DENDRO_573 + DENDRO_641) +
              DENDRO_160 * (DENDRO_674 + DENDRO_678) +
              DENDRO_160 * (-DENDRO_102 * DENDRO_94 + DENDRO_677) +
              DENDRO_160 * (DENDRO_191 * DENDRO_97 + DENDRO_680) +
              DENDRO_160 * (DENDRO_211 * DENDRO_55 + DENDRO_717) +
              DENDRO_234 * (DENDRO_458 + DENDRO_94 * grad_0_gt3[pp]) -
              DENDRO_261 * (DENDRO_354 + DENDRO_653) -
              DENDRO_267 * (DENDRO_421 + DENDRO_657) -
              DENDRO_267 * (DENDRO_655 + DENDRO_714) -
              DENDRO_267 * (DENDRO_423 + DENDRO_654 + DENDRO_715) -
              DENDRO_271 * (0.5 * DENDRO_456 + DENDRO_660) -
              DENDRO_271 * (DENDRO_454 + DENDRO_55 * DENDRO_94 + DENDRO_658) -
              DENDRO_681 * (DENDRO_350 + DENDRO_564 + DENDRO_595) + DENDRO_702 +
              DENDRO_718)) -
    DENDRO_20 *
        (DENDRO_713 +
         alpha[pp] *
             (DENDRO_117 * (DENDRO_405 + DENDRO_648) +
              DENDRO_117 * (DENDRO_405 + DENDRO_671) +
              DENDRO_117 * (DENDRO_409 + DENDRO_669) +
              DENDRO_117 * (DENDRO_473 + DENDRO_649) +
              DENDRO_117 * (DENDRO_566 + DENDRO_650) +
              DENDRO_117 * (DENDRO_646 + DENDRO_670) -
              DENDRO_154 * (DENDRO_451 + DENDRO_664) -
              DENDRO_154 * (DENDRO_505 + DENDRO_666) -
              DENDRO_154 * (DENDRO_508 + DENDRO_665) -
              DENDRO_154 * (DENDRO_514 + DENDRO_663) +
              DENDRO_160 * (DENDRO_675 + DENDRO_677) +
              DENDRO_160 * (DENDRO_678 + DENDRO_680) +
              DENDRO_160 * (DENDRO_229 * DENDRO_59 + DENDRO_674) +
              DENDRO_160 * (DENDRO_100 * DENDRO_94 + DENDRO_675 + DENDRO_676) +
              DENDRO_160 * (DENDRO_211 * DENDRO_56 + DENDRO_672 + DENDRO_673) +
              DENDRO_234 * (DENDRO_211 * grad_1_gt0[pp] + DENDRO_438) -
              DENDRO_261 * (DENDRO_352 + DENDRO_653) -
              DENDRO_261 * (DENDRO_173 * DENDRO_58 + DENDRO_596) -
              DENDRO_267 * (1.0 * DENDRO_428 + DENDRO_654) -
              DENDRO_267 * (0.5 * DENDRO_431 + DENDRO_657) -
              DENDRO_267 * (DENDRO_100 * DENDRO_211 + DENDRO_426 + DENDRO_655) -
              DENDRO_271 * (DENDRO_452 + DENDRO_660) -
              DENDRO_271 * (0.5 * DENDRO_453 + DENDRO_658) -
              DENDRO_271 * (DENDRO_100 * DENDRO_62 + DENDRO_287 + DENDRO_659) -
              DENDRO_651 * (DENDRO_113 + DENDRO_460 + DENDRO_461) -
              DENDRO_651 * (DENDRO_464 + DENDRO_465 + DENDRO_620) -
              DENDRO_681 * (DENDRO_378 + DENDRO_494 + DENDRO_577) +
              DENDRO_702)) -
    DENDRO_22 *
        (DENDRO_638 +
         alpha[pp] *
             (DENDRO_117 * (DENDRO_602 + DENDRO_645) +
              DENDRO_117 * (DENDRO_604 + DENDRO_608) +
              DENDRO_117 * (DENDRO_609 - DENDRO_643) +
              DENDRO_117 * (-DENDRO_183 * DENDRO_211 + DENDRO_644) +
              DENDRO_117 * (DENDRO_190 * DENDRO_226 + DENDRO_601 + DENDRO_645) -
              DENDRO_154 * (DENDRO_571 + DENDRO_594) -
              DENDRO_154 * (DENDRO_597 + DENDRO_642) -
              DENDRO_154 * (DENDRO_600 + DENDRO_642) +
              DENDRO_160 * (DENDRO_404 + DENDRO_612) +
              DENDRO_160 * (DENDRO_407 + DENDRO_619) +
              DENDRO_160 * (DENDRO_407 + DENDRO_650) +
              DENDRO_160 * (DENDRO_471 + DENDRO_648) +
              DENDRO_160 * (DENDRO_472 + DENDRO_649) +
              DENDRO_160 * (DENDRO_615 + DENDRO_646) -
              DENDRO_261 * (DENDRO_344 + DENDRO_587) -
              DENDRO_261 * (0.5 * DENDRO_346 + DENDRO_586) -
              DENDRO_267 * (DENDRO_590 + DENDRO_640) -
              DENDRO_267 * (DENDRO_591 + DENDRO_639) -
              DENDRO_267 * (DENDRO_190 * DENDRO_211 + DENDRO_419 + DENDRO_588) -
              DENDRO_271 * (DENDRO_451 + DENDRO_593) -
              DENDRO_271 * (DENDRO_118 * DENDRO_125 + DENDRO_641) +
              DENDRO_364 * (DENDRO_211 * grad_1_gt5[pp] + DENDRO_439) +
              DENDRO_629 - DENDRO_651 * (DENDRO_565 + DENDRO_595) +
              DENDRO_652)) -
    DENDRO_22 *
        (DENDRO_638 +
         alpha[pp] *
             (DENDRO_117 * (DENDRO_606 + DENDRO_607) +
              DENDRO_117 * (DENDRO_608 + DENDRO_609) +
              DENDRO_117 * (-DENDRO_102 * DENDRO_185 + DENDRO_604) +
              DENDRO_117 * (DENDRO_180 * DENDRO_211 + DENDRO_606) +
              DENDRO_117 * (DENDRO_191 * DENDRO_226 + DENDRO_602) -
              DENDRO_154 * (DENDRO_351 + DENDRO_594) -
              DENDRO_154 * (DENDRO_381 + DENDRO_497) -
              DENDRO_154 * (DENDRO_381 + DENDRO_598) -
              DENDRO_154 * (DENDRO_491 + DENDRO_596) -
              DENDRO_154 * (DENDRO_493 + DENDRO_597) -
              DENDRO_154 * (DENDRO_504 + DENDRO_600) +
              DENDRO_159 * (DENDRO_518 + DENDRO_584) +
              DENDRO_159 * (DENDRO_576 + DENDRO_583) +
              DENDRO_160 * (DENDRO_472 + DENDRO_614) +
              DENDRO_160 * (DENDRO_610 + DENDRO_612) +
              DENDRO_160 * (DENDRO_614 + DENDRO_615) +
              DENDRO_160 * (DENDRO_616 + DENDRO_619) -
              DENDRO_261 * (0.5 * DENDRO_363 + DENDRO_587) -
              DENDRO_261 * (1.0 * DENDRO_383 + DENDRO_585) -
              DENDRO_261 * (DENDRO_180 * DENDRO_226 + DENDRO_347 + DENDRO_586) -
              DENDRO_267 * (DENDRO_413 + DENDRO_590) -
              DENDRO_267 * (0.5 * DENDRO_418 + DENDRO_588) -
              DENDRO_267 * (DENDRO_180 * DENDRO_229 + DENDRO_416 + DENDRO_591) -
              DENDRO_271 * (DENDRO_450 + DENDRO_593) -
              DENDRO_271 * (DENDRO_100 * DENDRO_105 + DENDRO_579) +
              DENDRO_364 * (DENDRO_226 * grad_2_gt3[pp] + DENDRO_373) -
              DENDRO_621 * (DENDRO_466 + DENDRO_620) + DENDRO_629)) +
    DENDRO_26 *
        (DENDRO_308 * (-DENDRO_390 + DENDRO_391 - DENDRO_392 - DENDRO_395) +
         DENDRO_312 * (DENDRO_192 + DENDRO_388 * DENDRO_90) +
         DENDRO_389 * (DENDRO_229 + DENDRO_388 * DENDRO_70) +
         alpha[pp] *
             (DENDRO_117 * (DENDRO_413 + DENDRO_414) +
              DENDRO_117 * (DENDRO_416 - DENDRO_417) +
              DENDRO_117 * (DENDRO_419 + 1.0 * DENDRO_420) -
              DENDRO_125 * DENDRO_435 - DENDRO_154 * (DENDRO_402 + DENDRO_404) -
              DENDRO_154 * (DENDRO_409 - 1.0 * DENDRO_410) -
              DENDRO_154 * (DENDRO_175 * DENDRO_406 + DENDRO_405) -
              DENDRO_154 * (DENDRO_220 * DENDRO_406 + DENDRO_407) +
              DENDRO_160 * (DENDRO_421 + DENDRO_422) +
              DENDRO_160 * (DENDRO_423 + DENDRO_424) +
              DENDRO_160 * (DENDRO_426 + 1.0 * DENDRO_427) -
              DENDRO_173 * DENDRO_433 - DENDRO_175 * DENDRO_434 -
              DENDRO_212 * DENDRO_374 * grad_1_gt3[pp] -
              DENDRO_220 * DENDRO_261 * (DENDRO_181 - DENDRO_353) +
              DENDRO_234 * (DENDRO_425 + DENDRO_427) +
              DENDRO_234 * (DENDRO_125 * grad_2_gt3[pp] + DENDRO_431) +
              DENDRO_234 * (DENDRO_175 * grad_0_gt3[pp] + DENDRO_428) +
              DENDRO_358 * grad_0_gt3[pp] + DENDRO_360 * grad_1_gt3[pp] +
              DENDRO_362 * grad_2_gt3[pp] +
              DENDRO_364 * (DENDRO_418 + DENDRO_420) +
              DENDRO_364 * (DENDRO_173 * grad_0_gt3[pp] + DENDRO_429) +
              DENDRO_364 * (DENDRO_220 * grad_2_gt3[pp] + DENDRO_430) +
              DENDRO_367 * DENDRO_394 + DENDRO_440 -
              DENDRO_99 *
                  (DENDRO_396 + DENDRO_398 * DENDRO_95 +
                   DENDRO_401 * DENDRO_98 +
                   DENDRO_96 * (DENDRO_209 + DENDRO_399 + DENDRO_400))) -
         4 * grad2_1_1_alpha[pp]) +
    DENDRO_37 *
        (DENDRO_563 +
         alpha[pp] *
             (DENDRO_116 * (DENDRO_350 + DENDRO_565) +
              DENDRO_117 * (DENDRO_489 + DENDRO_569) +
              DENDRO_117 * (DENDRO_492 + DENDRO_570) +
              DENDRO_117 * (DENDRO_499 + DENDRO_571) +
              DENDRO_117 * (DENDRO_501 + DENDRO_570) -
              DENDRO_154 * (DENDRO_480 + DENDRO_567) -
              DENDRO_154 * (-DENDRO_121 * DENDRO_94 + DENDRO_484) -
              DENDRO_154 * (DENDRO_226 * DENDRO_58 + DENDRO_567) +
              DENDRO_160 * (DENDRO_451 + DENDRO_511) +
              DENDRO_160 * (DENDRO_510 + DENDRO_572) +
              DENDRO_160 * (DENDRO_513 + DENDRO_574) +
              DENDRO_160 * (DENDRO_573 + DENDRO_574) -
              DENDRO_165 * (DENDRO_459 + DENDRO_94 * grad_0_gt5[pp]) -
              DENDRO_261 * (0.5 * DENDRO_340 + DENDRO_468) -
              DENDRO_261 * (DENDRO_185 * DENDRO_58 + DENDRO_338 + DENDRO_467) -
              DENDRO_267 * (DENDRO_125 * DENDRO_180 + DENDRO_566) -
              DENDRO_271 * (1.0 * DENDRO_457 + DENDRO_477) -
              DENDRO_271 * (DENDRO_447 + DENDRO_475 + DENDRO_58 * DENDRO_94) -
              DENDRO_519 * (DENDRO_403 + DENDRO_576) + DENDRO_549 +
              DENDRO_582)) +
    DENDRO_37 *
        (DENDRO_563 +
         alpha[pp] *
             (DENDRO_117 * (DENDRO_355 + DENDRO_491) +
              DENDRO_117 * (DENDRO_355 + DENDRO_499) +
              DENDRO_117 * (DENDRO_379 + DENDRO_489) +
              DENDRO_117 * (DENDRO_492 + DENDRO_493) +
              DENDRO_117 * (DENDRO_495 + DENDRO_497) +
              DENDRO_117 * (DENDRO_501 + DENDRO_504) -
              DENDRO_154 * (DENDRO_482 + DENDRO_484) -
              DENDRO_154 * (DENDRO_485 + DENDRO_486) -
              DENDRO_154 * (DENDRO_206 * DENDRO_56 + DENDRO_479) -
              DENDRO_154 * (DENDRO_118 * DENDRO_94 + DENDRO_482 + DENDRO_483) -
              DENDRO_154 * (DENDRO_226 * DENDRO_59 + DENDRO_480 + DENDRO_481) +
              DENDRO_159 * (DENDRO_460 + DENDRO_463) +
              DENDRO_159 * (DENDRO_464 + DENDRO_466) +
              DENDRO_160 * (DENDRO_450 + DENDRO_510) +
              DENDRO_160 * (DENDRO_505 + DENDRO_506) +
              DENDRO_160 * (DENDRO_508 + DENDRO_511) +
              DENDRO_160 * (DENDRO_513 + DENDRO_514 + DENDRO_515) -
              DENDRO_165 * (DENDRO_226 * grad_2_gt0[pp] + DENDRO_371) -
              DENDRO_261 * (1.0 * DENDRO_365 + DENDRO_467) -
              DENDRO_261 * (0.5 * DENDRO_385 + DENDRO_469) -
              DENDRO_261 * (DENDRO_118 * DENDRO_226 + DENDRO_341 + DENDRO_468) -
              DENDRO_267 * (DENDRO_436 + DENDRO_474) -
              DENDRO_267 * (DENDRO_173 * DENDRO_55 + DENDRO_471) -
              DENDRO_271 * (DENDRO_281 + DENDRO_476) -
              DENDRO_271 * (0.5 * DENDRO_446 + DENDRO_475) -
              DENDRO_271 * (DENDRO_118 * DENDRO_97 + DENDRO_445 + DENDRO_477) -
              DENDRO_519 * (DENDRO_516 + DENDRO_518) + DENDRO_549)) +
    DENDRO_40 *
        (DENDRO_309 * DENDRO_63 +
         DENDRO_311 * (-DENDRO_79 - DENDRO_81 + DENDRO_83 - DENDRO_91) +
         DENDRO_389 * (DENDRO_71 + DENDRO_97) +
         alpha[pp] *
             (-DENDRO_115 * DENDRO_194 * DENDRO_276 +
              DENDRO_117 * (-1.0 * DENDRO_106 + DENDRO_114) +
              DENDRO_117 * (-1.0 * DENDRO_126 + DENDRO_132) +
              DENDRO_117 * (DENDRO_130 * DENDRO_286 + DENDRO_449) +
              DENDRO_117 * (DENDRO_168 * DENDRO_55 + DENDRO_451) +
              DENDRO_117 * (DENDRO_175 * DENDRO_58 + DENDRO_450) -
              DENDRO_125 * DENDRO_268 -
              DENDRO_130 * DENDRO_261 * (DENDRO_119 - DENDRO_135) -
              DENDRO_154 * (DENDRO_447 + 1.0 * DENDRO_448) -
              DENDRO_154 * (-DENDRO_152 * DENDRO_97 + DENDRO_445) +
              DENDRO_160 * (DENDRO_454 + 1.0 * DENDRO_455) +
              DENDRO_160 * (DENDRO_137 * DENDRO_97 + DENDRO_452) -
              DENDRO_165 * (DENDRO_446 + DENDRO_448) -
              DENDRO_165 * (DENDRO_130 * grad_2_gt0[pp] + DENDRO_457) +
              DENDRO_196 * DENDRO_357 + DENDRO_215 * DENDRO_359 -
              DENDRO_231 * DENDRO_269 * DENDRO_303 + DENDRO_233 * DENDRO_361 +
              DENDRO_234 * (DENDRO_453 + DENDRO_455) +
              DENDRO_234 * (DENDRO_125 * grad_2_gt0[pp] + DENDRO_456) +
              DENDRO_259 * DENDRO_366 - DENDRO_273 * DENDRO_458 -
              DENDRO_275 * DENDRO_459 + DENDRO_304 -
              DENDRO_99 *
                  (DENDRO_443 * DENDRO_96 + DENDRO_444 * DENDRO_98 + DENDRO_93 +
                   DENDRO_95 * (DENDRO_441 + DENDRO_442 + DENDRO_82))) -
         4 * grad2_0_0_alpha[pp]) +
    DENDRO_46 *
        (DENDRO_307 * DENDRO_309 + DENDRO_312 * (DENDRO_185 + DENDRO_310) +
         DENDRO_319 * (-DENDRO_313 + DENDRO_314 - DENDRO_315 - DENDRO_318) +
         alpha[pp] *
             (DENDRO_117 * (DENDRO_344 + DENDRO_345) +
              DENDRO_117 * (DENDRO_347 + 1.0 * DENDRO_348) -
              DENDRO_154 * (DENDRO_335 + DENDRO_336) -
              DENDRO_154 * (DENDRO_341 + 1.0 * DENDRO_342) -
              DENDRO_154 * (DENDRO_185 * DENDRO_339 + DENDRO_338) +
              DENDRO_160 * (DENDRO_349 + DENDRO_351) +
              DENDRO_160 * (DENDRO_118 * DENDRO_220 + DENDRO_354) +
              DENDRO_160 * (DENDRO_130 * DENDRO_180 + DENDRO_352) +
              DENDRO_160 * (DENDRO_168 * DENDRO_356 + DENDRO_355) -
              DENDRO_165 * (DENDRO_340 + DENDRO_342) -
              DENDRO_165 * (DENDRO_168 * grad_0_gt5[pp] + DENDRO_365) -
              DENDRO_168 * DENDRO_370 - DENDRO_173 * DENDRO_369 -
              DENDRO_186 * DENDRO_303 * (DENDRO_118 + DENDRO_119) -
              DENDRO_227 * DENDRO_374 * grad_2_gt5[pp] -
              DENDRO_273 * DENDRO_373 + DENDRO_317 * DENDRO_367 +
              DENDRO_358 * grad_0_gt5[pp] + DENDRO_360 * grad_1_gt5[pp] +
              DENDRO_362 * grad_2_gt5[pp] +
              DENDRO_364 * (DENDRO_346 + DENDRO_348) +
              DENDRO_364 * (DENDRO_173 * grad_0_gt5[pp] + DENDRO_363) -
              DENDRO_371 * DENDRO_372 + DENDRO_387 -
              DENDRO_99 *
                  (DENDRO_320 + DENDRO_326 * DENDRO_95 +
                   DENDRO_328 * DENDRO_96 +
                   DENDRO_98 * (DENDRO_330 + DENDRO_331 + DENDRO_332))) -
         4 * grad2_2_2_alpha[pp]);
const double DENDRO_720 = DENDRO_32 * DENDRO_719;
const double DENDRO_721 = (1.0 / 12.0) * chi[pp];
const double DENDRO_722 = (1.0 / 3.0) * At1[pp];
const double DENDRO_723 = At4[pp] * DENDRO_22;
const double DENDRO_724 = At3[pp] * DENDRO_26;
const double DENDRO_725 = DENDRO_34 + DENDRO_723 - DENDRO_724;
const double DENDRO_726 = At3[pp] * DENDRO_20;
const double DENDRO_727 = At4[pp] * DENDRO_37;
const double DENDRO_728 = -At1[pp] * DENDRO_40 + DENDRO_726 - DENDRO_727;
const double DENDRO_729 = At4[pp] * DENDRO_46;
const double DENDRO_730 =
    -At1[pp] * DENDRO_37 + At3[pp] * DENDRO_22 - DENDRO_729;
const double DENDRO_731 = -DENDRO_66 + DENDRO_67 + DENDRO_68;
const double DENDRO_732 = 6.0 * grad_2_alpha[pp];
const double DENDRO_733 = 6.0 * grad_0_alpha[pp];
const double DENDRO_734 = 6.0 * grad_1_alpha[pp];
const double DENDRO_735 = DENDRO_103 * DENDRO_143;
const double DENDRO_736 = -DENDRO_170 + DENDRO_171 + DENDRO_172;
const double DENDRO_737 = DENDRO_736 * grad_2_gt0[pp];
const double DENDRO_738 = DENDRO_137 * DENDRO_143;
const double DENDRO_739 = DENDRO_737 + DENDRO_738;
const double DENDRO_740 = DENDRO_208 - DENDRO_209 + DENDRO_210;
const double DENDRO_741 = DENDRO_149 * DENDRO_411;
const double DENDRO_742 = -DENDRO_134 * DENDRO_415;
const double DENDRO_743 = DENDRO_143 * DENDRO_280;
const double DENDRO_744 =
    DENDRO_157 * DENDRO_286 + 0.25 * DENDRO_736 * grad_0_gt0[pp];
const double DENDRO_745 = DENDRO_145 + DENDRO_151;
const double DENDRO_746 = DENDRO_138 * DENDRO_617;
const double DENDRO_747 = -DENDRO_217 + DENDRO_218 + DENDRO_219;
const double DENDRO_748 = DENDRO_180 * DENDRO_75;
const double DENDRO_749 = DENDRO_150 * DENDRO_747 + DENDRO_748;
const double DENDRO_750 = DENDRO_137 * DENDRO_149;
const double DENDRO_751 = DENDRO_736 * grad_1_gt0[pp] + DENDRO_750;
const double DENDRO_752 = 1.0 * DENDRO_116;
const double DENDRO_753 = DENDRO_137 * DENDRO_747;
const double DENDRO_754 = 2.0 * DENDRO_195;
const double DENDRO_755 = 2.0 * DENDRO_214;
const double DENDRO_756 = 2.0 * DENDRO_232;
const double DENDRO_757 = DENDRO_258 * DENDRO_52;
const double DENDRO_758 = (1.0 / 3.0) * At2[pp];
const double DENDRO_759 = At5[pp] * DENDRO_22;
const double DENDRO_760 =
    At2[pp] * DENDRO_20 - At4[pp] * DENDRO_26 + DENDRO_759;
const double DENDRO_761 = At5[pp] * DENDRO_37;
const double DENDRO_762 =
    -At2[pp] * DENDRO_40 + At4[pp] * DENDRO_20 - DENDRO_761;
const double DENDRO_763 = -At5[pp] * DENDRO_46 + DENDRO_39 + DENDRO_723;
const double DENDRO_764 = DENDRO_109 * DENDRO_149;
const double DENDRO_765 = DENDRO_138 * grad_2_gt5[pp];
const double DENDRO_766 = DENDRO_223 - DENDRO_224 + DENDRO_225;
const double DENDRO_767 = DENDRO_143 * DENDRO_337;
const double DENDRO_768 = DENDRO_178 - DENDRO_179 + DENDRO_184;
const double DENDRO_769 = DENDRO_138 * DENDRO_353;
const double DENDRO_770 = 0.25 * DENDRO_134 * grad_2_gt5[pp];
const double DENDRO_771 = DENDRO_135 * DENDRO_138;
const double DENDRO_772 = DENDRO_177 * DENDRO_75;
const double DENDRO_773 = DENDRO_149 * DENDRO_161;
const double DENDRO_774 = DENDRO_149 * DENDRO_337;
const double DENDRO_775 = -DENDRO_770;
const double DENDRO_776 = DENDRO_149 * grad_0_gt5[pp];
const double DENDRO_777 = 2 * At4[pp];
const double DENDRO_778 = At3[pp] * DENDRO_42;
const double DENDRO_779 = DENDRO_32 * DENDRO_777;
const double DENDRO_780 = 0.5 * DENDRO_388;
const double DENDRO_781 = DENDRO_32 * DENDRO_92;
const double DENDRO_782 = DENDRO_740 * grad_2_gt3[pp];
const double DENDRO_783 = DENDRO_134 * DENDRO_353;
const double DENDRO_784 = 0.25 * DENDRO_750;
const double DENDRO_785 = 1.0 * DENDRO_736;
const double DENDRO_786 = DENDRO_740 * grad_0_gt3[pp];
const double DENDRO_787 = (1.0 / 3.0) * At4[pp];
const double DENDRO_788 = DENDRO_747 * grad_2_gt5[pp];
const double DENDRO_789 = DENDRO_135 * DENDRO_747;
const double DENDRO_790 = -DENDRO_353 * DENDRO_747 + DENDRO_601;
const double DENDRO_791 = DENDRO_618 - DENDRO_783;
const double DENDRO_792 = DENDRO_766 * grad_1_gt5[pp];
const double DENDRO_793 = DENDRO_766 * grad_0_gt5[pp];
const double DENDRO_794 =
    DENDRO_321 * grad_2_chi[pp] + DENDRO_324 * grad_0_chi[pp] + DENDRO_86;
const double DENDRO_795 = 0.5 * DENDRO_306;
const double DENDRO_796 = DENDRO_32 * grad_0_alpha[pp];
const double DENDRO_797 = DENDRO_327 * grad_1_chi[pp] + DENDRO_49;
const double DENDRO_798 = DENDRO_32 * grad_1_alpha[pp];
const double DENDRO_799 =
    DENDRO_321 * grad_0_chi[pp] + DENDRO_329 * grad_2_chi[pp] + DENDRO_66;
const double DENDRO_800 = 0.5 * DENDRO_799;
const double DENDRO_801 = DENDRO_32 * grad_2_alpha[pp];
const double DENDRO_802 = 0.5 * DENDRO_797;
const double DENDRO_803 = 0.5 * grad_1_alpha[pp];
const double DENDRO_804 = 0.5 * grad_2_alpha[pp];
const double DENDRO_805 = DENDRO_244 * DENDRO_32;
const double DENDRO_806 = 0.5 * grad_0_alpha[pp];
const double DENDRO_807 = DENDRO_242 * DENDRO_32;
const double DENDRO_808 = DENDRO_239 * DENDRO_32;
const double DENDRO_809 = pow(DENDRO_20, 2);
const double DENDRO_810 = pow(DENDRO_37, 2);
const double DENDRO_811 = 2 * DENDRO_40;
const double DENDRO_812 = At0[pp] * pow(DENDRO_40, 2) + At3[pp] * DENDRO_809 +
                          At5[pp] * DENDRO_810 - DENDRO_239 * DENDRO_727 -
                          DENDRO_34 * DENDRO_811 + DENDRO_38 * DENDRO_811;
const double DENDRO_813 = 3 * DENDRO_115;
const double DENDRO_814 = pow(DENDRO_22, 2);
const double DENDRO_815 = 2 * DENDRO_26;
const double DENDRO_816 = At0[pp] * DENDRO_809 + At3[pp] * pow(DENDRO_26, 2) +
                          At5[pp] * DENDRO_814 + DENDRO_23 * DENDRO_239 -
                          DENDRO_34 * DENDRO_815 - DENDRO_723 * DENDRO_815;
const double DENDRO_817 = 2 * DENDRO_46;
const double DENDRO_818 = At0[pp] * DENDRO_810 + At3[pp] * DENDRO_814 +
                          At5[pp] * pow(DENDRO_46, 2) - DENDRO_242 * DENDRO_44 +
                          DENDRO_38 * DENDRO_817 - DENDRO_723 * DENDRO_817;
const double DENDRO_819 =
    At2[pp] * DENDRO_810 - DENDRO_20 * DENDRO_729 + DENDRO_22 * DENDRO_726 -
    DENDRO_22 * DENDRO_727 - DENDRO_34 * DENDRO_37 + DENDRO_37 * DENDRO_41 -
    DENDRO_40 * DENDRO_44 + DENDRO_40 * DENDRO_47 + DENDRO_46 * DENDRO_761;
const double DENDRO_820 = 6 * DENDRO_115;
const double DENDRO_821 = At1[pp] * DENDRO_809;
const double DENDRO_822 = DENDRO_20 * DENDRO_723;
const double DENDRO_823 = DENDRO_20 * DENDRO_38;
const double DENDRO_824 = DENDRO_26 * DENDRO_727;
const double DENDRO_825 = DENDRO_22 * DENDRO_761;
const double DENDRO_826 = DENDRO_20 * DENDRO_41;
const double DENDRO_827 = DENDRO_27 * DENDRO_40;
const double DENDRO_828 = DENDRO_23 * DENDRO_40;
const double DENDRO_829 = DENDRO_26 * DENDRO_726;
const double DENDRO_830 = DENDRO_821 + DENDRO_822 - DENDRO_823 + DENDRO_824 -
                          DENDRO_825 - DENDRO_826 + DENDRO_827 - DENDRO_828 -
                          DENDRO_829;
const double DENDRO_831 = At4[pp] * DENDRO_814;
const double DENDRO_832 = DENDRO_22 * DENDRO_34;
const double DENDRO_833 = DENDRO_21 * DENDRO_37;
const double DENDRO_834 = DENDRO_27 * DENDRO_37;
const double DENDRO_835 = DENDRO_22 * DENDRO_38;
const double DENDRO_836 = DENDRO_20 * DENDRO_47;
const double DENDRO_837 = DENDRO_22 * DENDRO_724;
const double DENDRO_838 = DENDRO_26 * DENDRO_729;
const double DENDRO_839 = DENDRO_46 * DENDRO_759;
const double DENDRO_840 = DENDRO_831 + DENDRO_832 - DENDRO_833 + DENDRO_834 -
                          DENDRO_835 - DENDRO_836 - DENDRO_837 + DENDRO_838 -
                          DENDRO_839;
const double DENDRO_841 = (1.0 / 3.0) * alpha[pp];
const double DENDRO_842 = 2 * DENDRO_115;
const double DENDRO_843 = DENDRO_842 * grad_0_alpha[pp];
const double DENDRO_844 = DENDRO_842 * grad_1_alpha[pp];
const double DENDRO_845 = DENDRO_842 * grad_2_alpha[pp];
const double DENDRO_846 = (7.0 / 3.0) * DENDRO_282;
const double DENDRO_847 = 4 * grad_0_K[pp];
const double DENDRO_848 = 9 * DENDRO_52;
const double DENDRO_849 = DENDRO_848 * DENDRO_95;
const double DENDRO_850 = DENDRO_32 * DENDRO_841;
const double DENDRO_851 = 4 * grad_2_K[pp];
const double DENDRO_852 = DENDRO_848 * DENDRO_98;
const double DENDRO_853 = 4 * grad_1_K[pp];
const double DENDRO_854 = DENDRO_848 * DENDRO_96;
const double DENDRO_855 = (1.0 / 3.0) * DENDRO_288;
const double DENDRO_856 = (1.0 / 3.0) * DENDRO_282;
const double DENDRO_857 = DENDRO_20 * DENDRO_32;
const double DENDRO_858 = (1.0 / 3.0) * DENDRO_857;
const double DENDRO_859 = (2.0 / 3.0) * DENDRO_16;
const double DENDRO_860 = (7.0 / 3.0) * DENDRO_857;
const double DENDRO_861 = pow(DENDRO_31, -3);
const double DENDRO_862 = DENDRO_0 * DENDRO_861;
const double DENDRO_863 = DENDRO_816 * DENDRO_862;
const double DENDRO_864 = DENDRO_818 * DENDRO_862;
const double DENDRO_865 = DENDRO_812 * DENDRO_862;
const double DENDRO_866 = 2.0 * DENDRO_861 * alpha[pp];
const double DENDRO_867 = DENDRO_830 * DENDRO_866;
const double DENDRO_868 = DENDRO_840 * DENDRO_866;
const double DENDRO_869 = DENDRO_819 * DENDRO_866;
const double DENDRO_870 = beta0[pp] * grad_0_Gt0[pp] +
                          beta1[pp] * grad_1_Gt0[pp] +
                          beta2[pp] * grad_2_Gt0[pp];
const double DENDRO_871 =
    DENDRO_168 * DENDRO_869 + DENDRO_173 * DENDRO_868 +
    DENDRO_175 * DENDRO_867 + DENDRO_185 * DENDRO_864 +
    DENDRO_192 * DENDRO_863 - 4.0 / 3.0 * DENDRO_288 * grad2_0_0_beta0[pp] -
    DENDRO_290 * grad2_1_1_beta0[pp] - DENDRO_292 * grad2_2_2_beta0[pp] +
    DENDRO_357 * DENDRO_859 - DENDRO_357 * grad_0_beta0[pp] -
    DENDRO_359 * grad_1_beta0[pp] - DENDRO_361 * grad_2_beta0[pp] +
    DENDRO_805 * grad2_1_2_beta0[pp] - DENDRO_812 * DENDRO_843 -
    DENDRO_819 * DENDRO_845 - DENDRO_830 * DENDRO_844 -
    DENDRO_846 * grad2_0_2_beta0[pp] -
    DENDRO_850 * (DENDRO_20 * DENDRO_853 + DENDRO_830 * DENDRO_854) -
    DENDRO_850 * (-DENDRO_37 * DENDRO_851 + DENDRO_819 * DENDRO_852) -
    DENDRO_850 * (-DENDRO_40 * DENDRO_847 + DENDRO_812 * DENDRO_849) -
    DENDRO_855 * grad2_0_1_beta1[pp] - DENDRO_855 * grad2_0_2_beta2[pp] -
    DENDRO_856 * grad2_1_2_beta1[pp] - DENDRO_856 * grad2_2_2_beta2[pp] +
    DENDRO_858 * grad2_1_1_beta1[pp] + DENDRO_858 * grad2_1_2_beta2[pp] +
    DENDRO_860 * grad2_0_1_beta0[pp] + DENDRO_865 * DENDRO_94 + DENDRO_870;
const double DENDRO_872 = -DENDRO_821 - DENDRO_822 + DENDRO_823 - DENDRO_824 +
                          DENDRO_825 + DENDRO_826 - DENDRO_827 + DENDRO_828 +
                          DENDRO_829;
const double DENDRO_873 = -DENDRO_831 - DENDRO_832 + DENDRO_833 - DENDRO_834 +
                          DENDRO_835 + DENDRO_836 + DENDRO_837 - DENDRO_838 +
                          DENDRO_839;
const double DENDRO_874 = DENDRO_20 * DENDRO_847;
const double DENDRO_875 = DENDRO_22 * DENDRO_851;
const double DENDRO_876 = DENDRO_26 * DENDRO_853;
const double DENDRO_877 = DENDRO_816 * DENDRO_854;
const double DENDRO_878 = (1.0 / 3.0) * DENDRO_290;
const double DENDRO_879 = DENDRO_22 * DENDRO_32;
const double DENDRO_880 = (1.0 / 3.0) * grad2_0_2_beta0[pp];
const double DENDRO_881 = (1.0 / 3.0) * DENDRO_879;
const double DENDRO_882 = (7.0 / 3.0) * DENDRO_879;
const double DENDRO_883 = beta0[pp] * grad_0_Gt1[pp] +
                          beta1[pp] * grad_1_Gt1[pp] +
                          beta2[pp] * grad_2_Gt1[pp];
const double DENDRO_884 =
    DENDRO_105 * DENDRO_869 + DENDRO_206 * DENDRO_864 -
    DENDRO_288 * grad2_0_0_beta1[pp] -
    4.0 / 3.0 * DENDRO_290 * grad2_1_1_beta1[pp] -
    DENDRO_292 * grad2_2_2_beta1[pp] + DENDRO_62 * DENDRO_865 -
    DENDRO_807 * grad2_0_2_beta1[pp] - DENDRO_816 * DENDRO_844 +
    DENDRO_858 * grad2_0_0_beta0[pp] + DENDRO_858 * grad2_0_2_beta2[pp] +
    DENDRO_860 * grad2_0_1_beta1[pp] - DENDRO_878 * grad2_0_1_beta0[pp] -
    DENDRO_878 * grad2_1_2_beta2[pp] + DENDRO_879 * DENDRO_880 +
    DENDRO_881 * grad2_2_2_beta2[pp] + DENDRO_882 * grad2_1_2_beta1[pp] +
    DENDRO_883;
const double DENDRO_885 = beta0[pp] * grad_0_Gt2[pp] +
                          beta1[pp] * grad_1_Gt2[pp] +
                          beta2[pp] * grad_2_Gt2[pp];
const double DENDRO_886 =
    DENDRO_125 * DENDRO_867 + DENDRO_130 * DENDRO_869 +
    DENDRO_220 * DENDRO_868 + DENDRO_226 * DENDRO_864 +
    DENDRO_229 * DENDRO_863 - DENDRO_288 * grad2_0_0_beta2[pp] -
    DENDRO_290 * grad2_1_1_beta2[pp] - DENDRO_292 * DENDRO_880 -
    1.0 / 3.0 * DENDRO_292 * grad2_1_2_beta1[pp] -
    4.0 / 3.0 * DENDRO_292 * grad2_2_2_beta2[pp] -
    DENDRO_357 * grad_0_beta2[pp] - DENDRO_359 * grad_1_beta2[pp] +
    DENDRO_361 * DENDRO_859 - DENDRO_361 * grad_2_beta2[pp] +
    DENDRO_808 * grad2_0_1_beta2[pp] - DENDRO_818 * DENDRO_845 -
    DENDRO_819 * DENDRO_843 - DENDRO_840 * DENDRO_844 -
    DENDRO_846 * grad2_0_2_beta2[pp] -
    DENDRO_850 * (DENDRO_22 * DENDRO_853 + DENDRO_840 * DENDRO_854) -
    DENDRO_850 * (-DENDRO_37 * DENDRO_847 + DENDRO_819 * DENDRO_849) -
    DENDRO_850 * (-DENDRO_46 * DENDRO_851 + DENDRO_818 * DENDRO_852) -
    DENDRO_856 * grad2_0_0_beta0[pp] - DENDRO_856 * grad2_0_1_beta1[pp] +
    DENDRO_865 * DENDRO_97 + DENDRO_881 * grad2_0_1_beta0[pp] +
    DENDRO_881 * grad2_1_1_beta1[pp] + DENDRO_882 * grad2_1_2_beta2[pp] +
    DENDRO_885;

// Dendro: printing variables
//--
a_rhs[pp]    = -DENDRO_0 * K[pp] + lambda[0] * (beta0[pp] * grad_0_alpha[pp] +
                                             beta1[pp] * grad_1_alpha[pp] +
                                             beta2[pp] * grad_2_alpha[pp]);
//--
b_rhs0[pp]   = B0[pp] * DENDRO_1 + lambda[1] * (beta0[pp] * grad_0_beta0[pp] +
                                              beta1[pp] * grad_1_beta0[pp] +
                                              beta2[pp] * grad_2_beta0[pp]);
//--
b_rhs1[pp]   = B1[pp] * DENDRO_1 + lambda[1] * (beta0[pp] * grad_0_beta1[pp] +
                                              beta1[pp] * grad_1_beta1[pp] +
                                              beta2[pp] * grad_2_beta1[pp]);
//--
b_rhs2[pp]   = B2[pp] * DENDRO_1 + lambda[1] * (beta0[pp] * grad_0_beta2[pp] +
                                              beta1[pp] * grad_1_beta2[pp] +
                                              beta2[pp] * grad_2_beta2[pp]);
//--
gt_rhs00[pp] = -At0[pp] * DENDRO_0 + DENDRO_2 * gt0[pp] +
               DENDRO_3 * grad_0_beta1[pp] + DENDRO_4 * grad_0_beta2[pp] -
               DENDRO_5 * grad_1_beta1[pp] - DENDRO_5 * grad_2_beta2[pp] +
               beta0[pp] * grad_0_gt0[pp] + beta1[pp] * grad_1_gt0[pp] +
               beta2[pp] * grad_2_gt0[pp];
//--
gt_rhs01[pp] = -At1[pp] * DENDRO_0 + DENDRO_6 * grad_0_beta0[pp] +
               DENDRO_6 * grad_1_beta1[pp] - DENDRO_7 * gt1[pp] +
               beta0[pp] * grad_0_gt1[pp] + beta1[pp] * grad_1_gt1[pp] +
               beta2[pp] * grad_2_gt1[pp] + grad_0_beta1[pp] * gt3[pp] +
               grad_0_beta2[pp] * gt4[pp] + grad_1_beta0[pp] * gt0[pp] +
               grad_1_beta2[pp] * gt2[pp];
//--
gt_rhs02[pp] = -At2[pp] * DENDRO_0 + DENDRO_8 * grad_0_beta0[pp] +
               DENDRO_8 * grad_2_beta2[pp] - DENDRO_9 * gt2[pp] +
               beta0[pp] * grad_0_gt2[pp] + beta1[pp] * grad_1_gt2[pp] +
               beta2[pp] * grad_2_gt2[pp] + grad_0_beta1[pp] * gt4[pp] +
               grad_0_beta2[pp] * gt5[pp] + grad_2_beta0[pp] * gt0[pp] +
               grad_2_beta1[pp] * gt1[pp];
//--
gt_rhs11[pp] = -At3[pp] * DENDRO_0 - DENDRO_10 * gt3[pp] + DENDRO_11 * gt3[pp] +
               DENDRO_12 * grad_1_beta2[pp] + DENDRO_3 * grad_1_beta0[pp] -
               DENDRO_7 * gt3[pp] + beta0[pp] * grad_0_gt3[pp] +
               beta1[pp] * grad_1_gt3[pp] + beta2[pp] * grad_2_gt3[pp];
//--
gt_rhs12[pp] = -At4[pp] * DENDRO_0 - DENDRO_10 * gt4[pp] +
               DENDRO_13 * grad_1_beta1[pp] + DENDRO_13 * grad_2_beta2[pp] +
               beta0[pp] * grad_0_gt4[pp] + beta1[pp] * grad_1_gt4[pp] +
               beta2[pp] * grad_2_gt4[pp] + grad_1_beta0[pp] * gt2[pp] +
               grad_1_beta2[pp] * gt5[pp] + grad_2_beta0[pp] * gt1[pp] +
               grad_2_beta1[pp] * gt3[pp];
//--
gt_rhs22[pp] = -At5[pp] * DENDRO_0 - DENDRO_10 * gt5[pp] +
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
    DENDRO_721 *
        (DENDRO_305 *
             (-DENDRO_117 * (DENDRO_106 - DENDRO_114) -
              DENDRO_117 * (DENDRO_126 + DENDRO_133) -
              DENDRO_117 * (DENDRO_136 + DENDRO_137 * DENDRO_139) -
              DENDRO_117 * (DENDRO_145 + DENDRO_149 * DENDRO_58) -
              DENDRO_117 * (DENDRO_143 * DENDRO_55 + DENDRO_151) +
              DENDRO_134 * DENDRO_268 -
              DENDRO_138 * DENDRO_261 * (DENDRO_120 + DENDRO_135) +
              DENDRO_154 * (DENDRO_163 + 1.0 * DENDRO_164) -
              DENDRO_154 * (-DENDRO_139 * DENDRO_59 + DENDRO_152 * DENDRO_75) +
              DENDRO_157 * DENDRO_270 * DENDRO_276 -
              DENDRO_160 * (DENDRO_156 + 1.0 * DENDRO_158) -
              DENDRO_160 *
                  (DENDRO_134 * DENDRO_161 + 1.0 * DENDRO_137 * DENDRO_75) +
              DENDRO_165 * (DENDRO_162 + DENDRO_164) +
              DENDRO_165 * (DENDRO_166 + DENDRO_167) - DENDRO_195 * DENDRO_196 -
              DENDRO_214 * DENDRO_215 - DENDRO_232 * DENDRO_233 -
              DENDRO_234 * (DENDRO_155 + DENDRO_158) -
              DENDRO_234 * (DENDRO_235 + DENDRO_236) - DENDRO_258 * DENDRO_259 +
              DENDRO_269 * DENDRO_271 * DENDRO_75 + DENDRO_272 * DENDRO_273 +
              DENDRO_274 * DENDRO_275 + DENDRO_304 -
              DENDRO_99 * (DENDRO_62 * DENDRO_96 + DENDRO_93 +
                           DENDRO_94 * DENDRO_95 + DENDRO_97 * DENDRO_98)) +
         DENDRO_63 * DENDRO_65 + DENDRO_720 * gt0[pp] -
         DENDRO_77 * (-DENDRO_71 + DENDRO_75) -
         DENDRO_92 * (DENDRO_79 + DENDRO_81 - DENDRO_83 + DENDRO_91) -
         12 * grad2_0_0_alpha[pp]) -
    alpha[pp] *
        (-At0[pp] * K[pp] + DENDRO_33 * (DENDRO_21 + DENDRO_23 - DENDRO_27) +
         DENDRO_43 * (DENDRO_34 + DENDRO_39 - DENDRO_41) +
         DENDRO_48 * (-At0[pp] * DENDRO_37 + DENDRO_44 - DENDRO_47)) +
    beta0[pp] * grad_0_At0[pp] + beta1[pp] * grad_1_At0[pp] +
    beta2[pp] * grad_2_At0[pp];
//--
At_rhs01[pp] =
    At0[pp] * grad_1_beta0[pp] - At1[pp] * DENDRO_7 +
    At2[pp] * grad_1_beta2[pp] + At3[pp] * grad_0_beta1[pp] +
    At4[pp] * grad_0_beta2[pp] +
    DENDRO_721 *
        (DENDRO_305 *
             (DENDRO_117 * (-DENDRO_286 * DENDRO_740 + DENDRO_716) -
              DENDRO_117 * (-DENDRO_616 + DENDRO_667 + DENDRO_668) +
              DENDRO_117 * (DENDRO_143 * DENDRO_411 - DENDRO_149 * DENDRO_617 +
                            DENDRO_647) +
              DENDRO_154 * (DENDRO_743 + DENDRO_744) +
              DENDRO_154 * (DENDRO_746 + DENDRO_749) +
              DENDRO_154 * (DENDRO_157 * DENDRO_406 + DENDRO_745) +
              DENDRO_154 * (DENDRO_133 + DENDRO_661 + DENDRO_662) +
              DENDRO_160 * (-DENDRO_55 * DENDRO_740 + DENDRO_717) +
              DENDRO_160 * (DENDRO_102 * DENDRO_157 - DENDRO_149 * DENDRO_280 +
                            DENDRO_676) -
              DENDRO_160 * (DENDRO_134 * DENDRO_498 + DENDRO_134 * DENDRO_512 +
                            DENDRO_190 * DENDRO_75) +
              DENDRO_160 * (-DENDRO_134 * DENDRO_617 - DENDRO_191 * DENDRO_75 +
                            DENDRO_679) -
              DENDRO_234 * (DENDRO_157 * grad_0_gt3[pp] + DENDRO_272) +
              DENDRO_260 * (DENDRO_735 + DENDRO_739) +
              DENDRO_261 * (-DENDRO_354 + DENDRO_500 + DENDRO_599) -
              DENDRO_267 * (DENDRO_102 * DENDRO_740 + DENDRO_714) -
              DENDRO_267 *
                  (-DENDRO_134 * DENDRO_380 + DENDRO_656 + DENDRO_742) -
              DENDRO_267 *
                  (-DENDRO_149 * DENDRO_284 + DENDRO_715 + DENDRO_741) +
              DENDRO_271 * (DENDRO_156 + DENDRO_157 * DENDRO_55 +
                            DENDRO_157 * DENDRO_56) +
              DENDRO_271 * (0.25 * DENDRO_235 + 0.5 * DENDRO_236 +
                            DENDRO_406 * DENDRO_75) +
              DENDRO_688 + DENDRO_689 + DENDRO_690 + DENDRO_691 + DENDRO_692 +
              DENDRO_693 + DENDRO_694 + DENDRO_695 + DENDRO_696 + DENDRO_697 +
              DENDRO_698 + DENDRO_699 + DENDRO_700 - DENDRO_701 * DENDRO_757 +
              DENDRO_718 -
              DENDRO_752 * (DENDRO_143 * grad_0_gt3[pp] + DENDRO_751) -
              DENDRO_752 * (DENDRO_134 * grad_1_gt5[pp] +
                            DENDRO_138 * grad_2_gt3[pp] + DENDRO_753) -
              DENDRO_754 * grad_0_gt1[pp] - DENDRO_755 * grad_1_gt1[pp] -
              DENDRO_756 * grad_2_gt1[pp] -
              DENDRO_99 * (DENDRO_112 * DENDRO_525 + DENDRO_125 * DENDRO_528 +
                           DENDRO_175 * DENDRO_523 + DENDRO_682)) +
         DENDRO_32 * DENDRO_732 * (DENDRO_125 - DENDRO_712 * DENDRO_731) +
         DENDRO_701 * DENDRO_719 -
         DENDRO_733 * (-DENDRO_703 + DENDRO_704 + DENDRO_705 + DENDRO_706) -
         DENDRO_734 * (DENDRO_707 + DENDRO_708 - DENDRO_709 - DENDRO_710) -
         12 * grad2_0_1_alpha[pp]) +
    DENDRO_722 * grad_0_beta0[pp] + DENDRO_722 * grad_1_beta1[pp] -
    alpha[pp] * (-At1[pp] * K[pp] + DENDRO_33 * DENDRO_725 +
                 DENDRO_43 * DENDRO_728 + DENDRO_48 * DENDRO_730) +
    beta0[pp] * grad_0_At1[pp] + beta1[pp] * grad_1_At1[pp] +
    beta2[pp] * grad_2_At1[pp];
//--
At_rhs02[pp] =
    At0[pp] * grad_2_beta0[pp] + At1[pp] * grad_2_beta1[pp] -
    At2[pp] * DENDRO_9 + At4[pp] * grad_0_beta1[pp] +
    At5[pp] * grad_0_beta2[pp] +
    DENDRO_721 *
        (DENDRO_305 *
             (-DENDRO_117 * (DENDRO_487 - DENDRO_488 + DENDRO_568) +
              DENDRO_117 * (-DENDRO_143 * DENDRO_498 - DENDRO_55 * DENDRO_768 +
                            DENDRO_774) -
              DENDRO_117 * (DENDRO_286 * DENDRO_766 + DENDRO_769 + DENDRO_770) +
              DENDRO_117 * (DENDRO_337 * DENDRO_747 - DENDRO_769 + DENDRO_775) -
              DENDRO_154 * (DENDRO_121 * DENDRO_139 - DENDRO_771 - DENDRO_772) -
              DENDRO_154 * (DENDRO_121 * DENDRO_157 - DENDRO_143 * DENDRO_161 -
                            DENDRO_60 * DENDRO_768) +
              DENDRO_154 * (DENDRO_58 * DENDRO_766 + DENDRO_771 + DENDRO_772) -
              DENDRO_160 * (DENDRO_136 + DENDRO_749) -
              DENDRO_160 * (DENDRO_744 + DENDRO_773) -
              DENDRO_160 * (DENDRO_157 * DENDRO_356 + DENDRO_745) -
              DENDRO_160 * (DENDRO_136 + DENDRO_138 * DENDRO_512 + DENDRO_748) +
              DENDRO_165 * (DENDRO_157 * grad_0_gt5[pp] + DENDRO_274) -
              DENDRO_261 * (DENDRO_121 * DENDRO_766 - 0.5 * DENDRO_765) -
              DENDRO_261 * (-DENDRO_135 * DENDRO_143 - DENDRO_58 * DENDRO_768 +
                            DENDRO_767) +
              DENDRO_266 * (DENDRO_751 + DENDRO_764) +
              DENDRO_267 * (DENDRO_134 * DENDRO_180 + 0.25 * DENDRO_753) +
              DENDRO_271 * (0.25 * DENDRO_166 + 1.0 * DENDRO_167) +
              DENDRO_271 * (DENDRO_157 * DENDRO_58 + DENDRO_157 * DENDRO_59 +
                            DENDRO_163) +
              DENDRO_530 + DENDRO_532 + DENDRO_534 + DENDRO_536 + DENDRO_538 +
              DENDRO_539 + DENDRO_541 + DENDRO_542 + DENDRO_543 + DENDRO_544 +
              DENDRO_545 + DENDRO_546 + DENDRO_547 - DENDRO_548 * DENDRO_757 +
              DENDRO_582 - DENDRO_752 * (DENDRO_739 + DENDRO_776) -
              DENDRO_754 * grad_0_gt2[pp] - DENDRO_755 * grad_1_gt2[pp] -
              DENDRO_756 * grad_2_gt2[pp] -
              DENDRO_99 * (DENDRO_105 * DENDRO_525 + DENDRO_130 * DENDRO_528 +
                           DENDRO_168 * DENDRO_523 + DENDRO_520)) +
         DENDRO_548 * DENDRO_719 + DENDRO_562 * DENDRO_734 -
         DENDRO_732 * (DENDRO_555 + DENDRO_556 - DENDRO_557 + DENDRO_558) -
         DENDRO_733 * (DENDRO_550 + DENDRO_551 - DENDRO_552 + DENDRO_553) -
         12 * grad2_0_2_alpha[pp]) +
    DENDRO_758 * grad_0_beta0[pp] + DENDRO_758 * grad_2_beta2[pp] -
    alpha[pp] * (-At2[pp] * K[pp] + DENDRO_33 * DENDRO_760 +
                 DENDRO_43 * DENDRO_762 + DENDRO_48 * DENDRO_763) +
    beta0[pp] * grad_0_At2[pp] + beta1[pp] * grad_1_At2[pp] +
    beta2[pp] * grad_2_At2[pp];
//--
At_rhs11[pp] =
    -At3[pp] * DENDRO_10 + At3[pp] * DENDRO_11 - At3[pp] * DENDRO_7 +
    DENDRO_17 * grad_1_beta0[pp] +
    DENDRO_721 *
        (DENDRO_305 *
             (-DENDRO_117 * (DENDRO_412 - 1.0 * DENDRO_414) +
              DENDRO_117 * (DENDRO_419 - 1.0 * DENDRO_782) -
              DENDRO_117 * (DENDRO_415 * DENDRO_747 + DENDRO_417) +
              DENDRO_134 * DENDRO_435 + DENDRO_149 * DENDRO_434 +
              DENDRO_154 * (-DENDRO_409 + DENDRO_410) +
              DENDRO_154 * (DENDRO_144 * DENDRO_736 + DENDRO_149 * DENDRO_406) +
              DENDRO_154 * (DENDRO_406 * DENDRO_747 + DENDRO_783) +
              DENDRO_154 * (DENDRO_56 * DENDRO_785 + DENDRO_784) +
              DENDRO_160 * (DENDRO_422 + DENDRO_742) +
              DENDRO_160 * (DENDRO_424 + DENDRO_741) +
              DENDRO_160 * (DENDRO_426 - 1.0 * DENDRO_786) +
              DENDRO_234 * (DENDRO_425 - DENDRO_786) +
              DENDRO_234 * (-DENDRO_134 * grad_2_gt3[pp] + DENDRO_431) +
              DENDRO_234 * (-DENDRO_149 * grad_0_gt3[pp] + DENDRO_428) -
              DENDRO_261 * DENDRO_747 * (DENDRO_182 + DENDRO_353) +
              6.0 * DENDRO_266 * DENDRO_740 * grad_1_gt3[pp] +
              DENDRO_364 * (DENDRO_418 - DENDRO_782) +
              DENDRO_364 * (DENDRO_429 - DENDRO_736 * grad_0_gt3[pp]) +
              DENDRO_364 * (DENDRO_430 - DENDRO_747 * grad_2_gt3[pp]) -
              DENDRO_394 * DENDRO_757 + DENDRO_433 * DENDRO_736 + DENDRO_440 -
              DENDRO_754 * grad_0_gt3[pp] - DENDRO_755 * grad_1_gt3[pp] -
              DENDRO_756 * grad_2_gt3[pp] -
              DENDRO_99 * (DENDRO_192 * DENDRO_95 + DENDRO_211 * DENDRO_96 +
                           DENDRO_229 * DENDRO_98 + DENDRO_396)) -
         DENDRO_64 * (DENDRO_390 - DENDRO_391 + DENDRO_392 + DENDRO_395) +
         DENDRO_720 * gt3[pp] +
         DENDRO_77 * (DENDRO_229 - DENDRO_731 * DENDRO_780) +
         DENDRO_781 *
             (DENDRO_192 - DENDRO_780 * (-DENDRO_86 + DENDRO_87 + DENDRO_88)) -
         12 * grad2_1_1_alpha[pp]) +
    DENDRO_777 * grad_1_beta2[pp] -
    alpha[pp] * (-At3[pp] * K[pp] + DENDRO_33 * DENDRO_728 +
                 DENDRO_725 * DENDRO_778 + DENDRO_730 * DENDRO_779) +
    beta0[pp] * grad_0_At3[pp] + beta1[pp] * grad_1_At3[pp] +
    beta2[pp] * grad_2_At3[pp];
//--
At_rhs12[pp] =
    At1[pp] * grad_2_beta0[pp] + At2[pp] * grad_1_beta0[pp] +
    At3[pp] * grad_2_beta1[pp] - At4[pp] * DENDRO_10 +
    At5[pp] * grad_1_beta2[pp] +
    DENDRO_721 *
        (DENDRO_305 *
             (DENDRO_117 * (DENDRO_183 * DENDRO_740 + DENDRO_644) +
              DENDRO_117 * (-DENDRO_190 * DENDRO_766 + DENDRO_790) +
              DENDRO_117 * (DENDRO_333 * DENDRO_747 + DENDRO_790) -
              DENDRO_117 * (DENDRO_100 * DENDRO_768 + DENDRO_498 * DENDRO_736 +
                            DENDRO_643) +
              DENDRO_117 * (-DENDRO_512 * DENDRO_736 + DENDRO_603 -
                            DENDRO_617 * DENDRO_736) +
              DENDRO_153 * (DENDRO_735 + DENDRO_737 + DENDRO_776) -
              DENDRO_154 * (DENDRO_139 * DENDRO_183 + DENDRO_775 - DENDRO_789) -
              DENDRO_154 * (-DENDRO_161 * DENDRO_736 - DENDRO_56 * DENDRO_768 +
                            DENDRO_774) +
              DENDRO_154 * (DENDRO_406 * DENDRO_766 + DENDRO_770 + DENDRO_789) +
              DENDRO_160 * (-DENDRO_138 * DENDRO_380 + DENDRO_791) +
              DENDRO_160 * (-DENDRO_356 * DENDRO_740 + DENDRO_646) +
              DENDRO_160 * (-DENDRO_617 * DENDRO_747 + DENDRO_791) +
              DENDRO_160 *
                  (-DENDRO_143 * DENDRO_284 + DENDRO_647 - 0.25 * DENDRO_764) +
              DENDRO_160 *
                  (-DENDRO_280 * DENDRO_736 + DENDRO_611 - DENDRO_784) +
              DENDRO_160 *
                  (-DENDRO_406 * DENDRO_740 + DENDRO_472 + DENDRO_613) -
              DENDRO_261 * (DENDRO_183 * DENDRO_766 - 0.5 * DENDRO_788) -
              DENDRO_261 * (-DENDRO_135 * DENDRO_736 + DENDRO_337 * DENDRO_736 -
                            DENDRO_406 * DENDRO_768) -
              DENDRO_267 * (-DENDRO_380 * DENDRO_747 + DENDRO_639) -
              DENDRO_267 * (-DENDRO_190 * DENDRO_740 - DENDRO_191 * DENDRO_740 +
                            DENDRO_419) -
              DENDRO_267 *
                  (-DENDRO_284 * DENDRO_736 + DENDRO_589 + DENDRO_640) +
              DENDRO_271 * (DENDRO_118 * DENDRO_134 + DENDRO_746) +
              DENDRO_271 * (DENDRO_151 + DENDRO_743 + DENDRO_773) +
              DENDRO_364 * (DENDRO_439 - DENDRO_740 * grad_1_gt5[pp]) -
              DENDRO_627 * DENDRO_757 + DENDRO_628 + DENDRO_652 -
              DENDRO_754 * grad_0_gt4[pp] - DENDRO_755 * grad_1_gt4[pp] -
              DENDRO_756 * grad_2_gt4[pp] -
              DENDRO_99 * (DENDRO_173 * DENDRO_523 + DENDRO_202 * DENDRO_525 +
                           DENDRO_220 * DENDRO_528 + DENDRO_622)) -
         DENDRO_32 * DENDRO_733 * (-DENDRO_637 + DENDRO_736) +
         DENDRO_627 * DENDRO_719 + DENDRO_631 * DENDRO_734 -
         DENDRO_732 * (-DENDRO_632 + DENDRO_633 + DENDRO_634 + DENDRO_635) -
         12 * grad2_1_2_alpha[pp]) +
    DENDRO_787 * grad_1_beta1[pp] + DENDRO_787 * grad_2_beta2[pp] -
    alpha[pp] * (-At4[pp] * K[pp] + DENDRO_33 * DENDRO_762 +
                 DENDRO_760 * DENDRO_778 + DENDRO_763 * DENDRO_779) +
    beta0[pp] * grad_0_At4[pp] + beta1[pp] * grad_1_At4[pp] +
    beta2[pp] * grad_2_At4[pp];
//--
At_rhs22[pp] =
    -At5[pp] * DENDRO_10 + At5[pp] * DENDRO_14 - At5[pp] * DENDRO_9 +
    DENDRO_18 * grad_2_beta0[pp] +
    DENDRO_721 *
        (DENDRO_305 *
             (-DENDRO_117 * (DENDRO_343 - 1.0 * DENDRO_345) -
              DENDRO_117 * (0.25 * DENDRO_788 + 1.0 * DENDRO_792) +
              DENDRO_138 * DENDRO_372 * grad_0_gt5[pp] +
              DENDRO_143 * DENDRO_370 +
              DENDRO_154 * (DENDRO_334 - 1.0 * DENDRO_336) +
              DENDRO_154 * (0.25 * DENDRO_765 + 1.0 * DENDRO_793) -
              DENDRO_154 * (-DENDRO_339 * DENDRO_768 + DENDRO_767) -
              DENDRO_160 * (DENDRO_118 * DENDRO_747 + DENDRO_769) -
              DENDRO_160 * (DENDRO_138 * DENDRO_180 + DENDRO_789) -
              DENDRO_160 * (DENDRO_143 * DENDRO_356 + DENDRO_150 * DENDRO_736) -
              DENDRO_160 * (DENDRO_59 * DENDRO_785 + 0.25 * DENDRO_738) +
              DENDRO_165 * (DENDRO_765 + DENDRO_793) +
              DENDRO_165 *
                  (DENDRO_143 * grad_0_gt5[pp] + DENDRO_768 * grad_2_gt0[pp]) +
              6.0 * DENDRO_260 * DENDRO_766 * grad_2_gt5[pp] -
              DENDRO_261 * DENDRO_768 * (DENDRO_120 + DENDRO_322) +
              DENDRO_273 * DENDRO_747 * grad_1_gt5[pp] -
              DENDRO_317 * DENDRO_757 - DENDRO_364 * (DENDRO_788 + DENDRO_792) -
              DENDRO_364 *
                  (DENDRO_137 * DENDRO_768 + DENDRO_736 * grad_0_gt5[pp]) +
              DENDRO_369 * DENDRO_736 + DENDRO_387 -
              DENDRO_754 * grad_0_gt5[pp] - DENDRO_755 * grad_1_gt5[pp] -
              DENDRO_756 * grad_2_gt5[pp] -
              DENDRO_99 * (DENDRO_185 * DENDRO_95 + DENDRO_206 * DENDRO_96 +
                           DENDRO_226 * DENDRO_98 + DENDRO_320)) +
         DENDRO_307 * DENDRO_65 + DENDRO_720 * gt5[pp] -
         DENDRO_76 * (DENDRO_313 - DENDRO_314 + DENDRO_315 + DENDRO_318) -
         DENDRO_781 * (-DENDRO_310 + DENDRO_768) - 12 * grad2_2_2_alpha[pp]) +
    DENDRO_777 * grad_2_beta1[pp] -
    alpha[pp] * (At5[pp] * DENDRO_42 * DENDRO_763 - At5[pp] * K[pp] +
                 DENDRO_48 * DENDRO_762 + DENDRO_760 * DENDRO_779) +
    beta0[pp] * grad_0_At5[pp] + beta1[pp] * grad_1_At5[pp] +
    beta2[pp] * grad_2_At5[pp];
//--
K_rhs[pp] =
    -DENDRO_288 * chi[pp] *
        (DENDRO_798 * (DENDRO_443 + DENDRO_53 * DENDRO_802) +
         DENDRO_801 * (DENDRO_444 + DENDRO_53 * DENDRO_800) -
         grad2_0_0_alpha[pp] +
         grad_0_alpha[pp] *
             (DENDRO_32 * DENDRO_441 + DENDRO_32 * DENDRO_442 -
              DENDRO_52 * (-0.5 * DENDRO_794 * DENDRO_85 + DENDRO_84) +
              DENDRO_83)) -
    DENDRO_290 * chi[pp] *
        (DENDRO_796 * (DENDRO_398 + DENDRO_780 * DENDRO_794) +
         DENDRO_801 * (DENDRO_401 + DENDRO_780 * DENDRO_799) -
         grad2_1_1_alpha[pp] +
         grad_1_alpha[pp] *
             (DENDRO_32 * DENDRO_399 + DENDRO_32 * DENDRO_400 + DENDRO_391 -
              DENDRO_52 * (DENDRO_393 - DENDRO_394 * DENDRO_802))) -
    DENDRO_292 * chi[pp] *
        (DENDRO_796 * (DENDRO_326 + DENDRO_794 * DENDRO_795) +
         DENDRO_798 * (DENDRO_328 + DENDRO_795 * DENDRO_797) -
         grad2_2_2_alpha[pp] +
         grad_2_alpha[pp] *
             (DENDRO_32 * DENDRO_330 + DENDRO_32 * DENDRO_331 +
              DENDRO_32 * DENDRO_332 -
              DENDRO_52 * (DENDRO_316 - DENDRO_317 * DENDRO_800))) +
    DENDRO_805 * chi[pp] *
        (0.5 * DENDRO_796 * (DENDRO_623 + DENDRO_636 * DENDRO_794) +
         DENDRO_803 * (DENDRO_32 * DENDRO_624 -
                       DENDRO_52 * (-DENDRO_627 * DENDRO_797 + grad_2_chi[pp]) +
                       DENDRO_630) +
         DENDRO_804 * (DENDRO_32 * DENDRO_625 + DENDRO_32 * DENDRO_626 -
                       DENDRO_52 * (-DENDRO_627 * DENDRO_799 + grad_1_chi[pp]) +
                       DENDRO_632) -
         grad2_1_2_alpha[pp]) -
    DENDRO_807 * chi[pp] *
        (0.5 * DENDRO_798 * (DENDRO_524 + DENDRO_561 * DENDRO_797) +
         DENDRO_804 * (DENDRO_32 * DENDRO_526 + DENDRO_32 * DENDRO_527 -
                       DENDRO_52 * (-DENDRO_548 * DENDRO_799 + grad_0_chi[pp]) +
                       DENDRO_557) +
         DENDRO_806 * (DENDRO_32 * DENDRO_521 + DENDRO_32 * DENDRO_522 -
                       DENDRO_52 * (-DENDRO_548 * DENDRO_794 + grad_2_chi[pp]) +
                       DENDRO_552) -
         grad2_0_2_alpha[pp]) +
    DENDRO_808 * chi[pp] *
        (0.5 * DENDRO_801 * (DENDRO_686 + DENDRO_712 * DENDRO_799) +
         DENDRO_803 * (DENDRO_32 * DENDRO_685 -
                       DENDRO_52 * (-DENDRO_701 * DENDRO_797 + grad_0_chi[pp]) +
                       DENDRO_711) +
         DENDRO_806 * (DENDRO_32 * DENDRO_683 + DENDRO_32 * DENDRO_684 -
                       DENDRO_52 * (-DENDRO_701 * DENDRO_794 + grad_1_chi[pp]) +
                       DENDRO_703) -
         grad2_0_1_alpha[pp]) +
    DENDRO_841 *
        (At0[pp] * DENDRO_812 * DENDRO_813 + At1[pp] * DENDRO_820 * DENDRO_830 +
         At2[pp] * DENDRO_819 * DENDRO_820 + At3[pp] * DENDRO_813 * DENDRO_816 +
         At4[pp] * DENDRO_820 * DENDRO_840 + At5[pp] * DENDRO_813 * DENDRO_818 +
         pow(K[pp], 2)) +
    beta0[pp] * grad_0_K[pp] + beta1[pp] * grad_1_K[pp] +
    beta2[pp] * grad_2_K[pp];
//--
Gt_rhs0[pp] = DENDRO_871;
//--
Gt_rhs1[pp] = -DENDRO_112 * DENDRO_866 * DENDRO_872 +
              DENDRO_195 * grad_0_beta1[pp] -
              DENDRO_202 * DENDRO_866 * DENDRO_873 - DENDRO_214 * DENDRO_859 +
              DENDRO_214 * grad_1_beta1[pp] + DENDRO_232 * grad_2_beta1[pp] -
              DENDRO_740 * DENDRO_863 + DENDRO_843 * DENDRO_872 +
              DENDRO_845 * DENDRO_873 + DENDRO_850 * (DENDRO_876 - DENDRO_877) -
              DENDRO_850 * (-DENDRO_849 * DENDRO_872 + DENDRO_874) -
              DENDRO_850 * (-DENDRO_852 * DENDRO_873 + DENDRO_875) + DENDRO_884;
//--
Gt_rhs2[pp] = DENDRO_886;
//--
B_rhs0[pp] =
    -B0[pp] * eta - DENDRO_870 * lambda[3] + DENDRO_871 +
    lambda[2] * (beta0[pp] * grad_0_B0[pp] + beta1[pp] * grad_1_B0[pp] +
                 beta2[pp] * grad_2_B0[pp]);
//--
B_rhs1[pp] =
    -B1[pp] * eta + DENDRO_112 * DENDRO_867 + DENDRO_202 * DENDRO_868 +
    DENDRO_211 * DENDRO_863 - DENDRO_357 * grad_0_beta1[pp] +
    DENDRO_359 * DENDRO_859 - DENDRO_359 * grad_1_beta1[pp] -
    DENDRO_361 * grad_2_beta1[pp] - DENDRO_830 * DENDRO_843 -
    DENDRO_840 * DENDRO_845 - DENDRO_850 * (-DENDRO_876 + DENDRO_877) -
    DENDRO_850 * (DENDRO_830 * DENDRO_849 + DENDRO_874) -
    DENDRO_850 * (DENDRO_840 * DENDRO_852 + DENDRO_875) -
    DENDRO_883 * lambda[3] + DENDRO_884 +
    lambda[2] * (beta0[pp] * grad_0_B1[pp] + beta1[pp] * grad_1_B1[pp] +
                 beta2[pp] * grad_2_B1[pp]);
//--
B_rhs2[pp] =
    -B2[pp] * eta - DENDRO_885 * lambda[3] + DENDRO_886 +
    lambda[2] * (beta0[pp] * grad_0_B2[pp] + beta1[pp] * grad_1_B2[pp] +
                 beta2[pp] * grad_2_B2[pp]);
// Dendro: reduced ops: 3952
// Dendro: }}}
