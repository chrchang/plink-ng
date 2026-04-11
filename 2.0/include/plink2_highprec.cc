// This library is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
// Christopher Chang.
//
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library.  If not, see <http://www.gnu.org/licenses/>.

#include "plink2_highprec.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
namespace plink2 {
#endif

// Portable "double-double" operations supporting high-accuracy log-likelihood
// calculations, based on a small subset of the QD library
// (https://github.com/BL-highprecision/QD ).  See LICENSE.QD for that
// library's BSD-3-Clause-LBNL license.

const dd_real _ddr_e = {{2.718281828459045091e+00, 1.445646891729250158e-16}};
const dd_real _ddr_log2 = {{6.931471805599452862e-01, 2.319046813846299558e-17}};
const dd_real _ddr_half_log_2pi = {{9.1893853320467278056e-01, -3.8782941580672414498e-17}};
const dd_real _ddr_12th = {{8.3333333333333328707e-02,  4.6259292692714853283e-18}};
const dd_real _ddr_1188th = {{8.4175084175084171397e-04,  3.6870174889237693563e-20}};
const dd_real _ddr_1260th = {{7.9365079365079365011e-04,  6.8838233173682821147e-22}};
const dd_real _ddr_neg_360th = {{-2.7777777777777778838e-03,  1.0601087908747154118e-19}};
const dd_real _ddr_neg_1680th = {{-5.9523809523809529179e-04,  5.3693821875472602376e-20}};

const double _ddr_eps = 4.93038065763132e-32;  // 2^-104

// These values were calculated with the quad-double side of the QD library as
// follows:
//   #include <qd/qd_real.h>
//   ...
//   qd_real log_sum = qd_real(0);
//   for (int32_t ii = 1; ii < 256; ii++) {
//     fprintf(stderr, "ln(%d!): ", ii);
//     log_sum += log(qd_real(ii));
//     log_sum.dump();
//   }
// and then postprocessed into the array declaration below.
const dd_real _ddr_ln_fact[_ddr_n_ln_fact] = {
  { 0.0000000000000000000e+00,  0.0000000000000000000e+00},
  { 0.0000000000000000000e+00,  0.0000000000000000000e+00},
  { 6.9314718055994528623e-01,  2.3190468138462995584e-17},
  { 1.7917594692280549573e+00,  4.3499798250963346992e-17},
  { 3.1780538303479457518e+00, -1.3216387039714196992e-16},
  { 4.7874917427820458116e+00,  1.8268155143874836655e-16},
  { 6.5792512120101012130e+00, -2.1790786016035090879e-16},
  { 8.5251613610654146669e+00, -3.6671660300633312458e-16},
  { 1.0604602902745250859e+01, -6.3021210597849109071e-16},
  { 1.2801827480081469091e+01,  5.2062957887166613404e-16},
  { 1.5104412573075515880e+01, -5.8462446316668400873e-16},
  { 1.7502307845873886549e+01, -7.0998288430900017421e-16},
  { 1.9987214495661884683e+01,  1.4661311288682236414e-15},
  { 2.2552163853123424531e+01, -1.6450514375919355420e-15},
  { 2.5191221182738679829e+01,  1.6710216640385303929e-15},
  { 2.7899271383840890337e+01,  1.2290202987493116658e-15},
  { 3.0671860106080671926e+01,  8.7769296145310094573e-16},
  { 3.3505073450136890756e+01, -1.8718492640011379814e-15},
  { 3.6395445208033052609e+01,  9.6751742592521715127e-16},
  { 3.9339884187199494647e+01, -6.1107769227967714765e-16},
  { 4.2335616460753485057e+01, -2.7806729241829511925e-17},
  { 4.5380138898476907627e+01,  3.9880537033726689082e-16},
  { 4.8471181351835227247e+01, -3.3670985639296026129e-15},
  { 5.1606675567764376922e+01, -3.3513402040623197993e-15},
  { 5.4784729398112318677e+01,  5.1329881419110192418e-16},
  { 5.8003605222980517908e+01,  2.0311680775630076815e-15},
  { 6.1261701761002001376e+01,  6.0851051617904650468e-16},
  { 6.4557538627006337606e+01, -6.5470111535469699323e-15},
  { 6.7889743137181540078e+01, -5.0951267256408069202e-15},
  { 7.1257038967168014665e+01, -5.6547469778977255102e-15},
  { 7.4658236348830158136e+01,  6.2499176982907561629e-15},
  { 7.8092223553315307072e+01,  3.5595186754960832134e-15},
  { 8.1557959456115042940e+01, -5.7614246931254325841e-15},
  { 8.5054467017581515620e+01,  1.7940650857583316722e-15},
  { 8.8580827542197681623e+01, -2.8196658134202104736e-15},
  { 9.2136175603687092917e+01, -4.3318308517998912416e-16},
  { 9.5719694542143201943e+01,  5.4199493102206283914e-16},
  { 9.9330612454787427623e+01, -6.9383351466841645262e-16},
  { 1.0296819861451380973e+02,  2.9688100510033878885e-15},
  { 1.0663176026064346047e+02, -1.3433085124319429993e-15},
  { 1.1032063971475739095e+02,  4.4812011344826036031e-15},
  { 1.1403421178146170689e+02, -3.6586846922455711238e-15},
  { 1.1777188139974506953e+02,  2.0091660912102236417e-15},
  { 1.2153308151543863858e+02, -4.6133769756837181586e-15},
  { 1.2531727114935689826e+02, -3.1380422260738893674e-15},
  { 1.2912393363912721611e+02, -1.2282659095377790667e-15},
  { 1.3295257503561632006e+02, -1.0182123580995801232e-14},
  { 1.3680272263732635452e+02,  1.3954035745823505920e-14},
  { 1.4067392364823425055e+02,  8.8490587327516230589e-15},
  { 1.4456574394634489522e+02, -9.2121271469428466604e-15},
  { 1.4847776695177302031e+02,  1.1757835515507762342e-14},
  { 1.5240959258449734648e+02,  1.1360070971878851950e-14},
  { 1.5636083630307879844e+02, -1.3243057336032416937e-14},
  { 1.6033112821663090131e+02,  5.7217979044139603522e-15},
  { 1.6432011226319517050e+02,  1.0912942276165645104e-14},
  { 1.6832744544842765322e+02, -8.8797938909247163255e-16},
  { 1.7235279713916278865e+02,  1.2910571080291391089e-14},
  { 1.7639584840699734514e+02,  6.5783262863108245512e-15},
  { 1.8045629141754378111e+02, -1.0056337354872399496e-14},
  { 1.8453382886144947861e+02,  1.1893610811486772714e-14},
  { 1.8862817342367159767e+02, -6.4876226164530563651e-15},
  { 1.9273904728784489748e+02,  4.9535392715078671495e-15},
  { 1.9686618167289000780e+02, -1.3811903140213112046e-14},
  { 2.0100931639928151640e+02,  1.0282768713994318190e-14},
  { 2.0516819948264119944e+02, -9.0236332835150033269e-16},
  { 2.0934258675253684601e+02, -1.0360038620726708595e-14},
  { 2.1353224149456326586e+02, -4.6687375155672468668e-15},
  { 2.1773693411395422004e+02,  7.2129167385410351478e-15},
  { 2.2195644181913033322e+02,  7.3499716563819027639e-16},
  { 2.2619054832372759734e+02, -4.0028941502701989522e-15},
  { 2.3043904356577695580e+02, -3.4806000957542813637e-15},
  { 2.3470172344281826327e+02,  4.4772149730333026359e-15},
  { 2.3897838956183431947e+02,  3.5882043155110511340e-15},
  { 2.4326884900298270509e+02,  9.0952962190379383621e-15},
  { 2.4757291409618687794e+02,  5.9952790996231564433e-15},
  { 2.5189040220972319162e+02,  2.7594986872193897642e-15},
  { 2.5622113555000953511e+02, -9.6529011360351126336e-15},
  { 2.6056494097186322279e+02, -1.3479781978823912067e-14},
  { 2.6492164979855277807e+02,  2.2976474929622466063e-14},
  { 2.6929109765101981111e+02,  1.1426270998958740688e-14},
  { 2.7367312428569368876e+02,  1.5386591972148984368e-14},
  { 2.7806757344036611812e+02,  2.4793702699450299690e-14},
  { 2.8247429268763039545e+02,  5.7877348379581861847e-16},
  { 2.8689313329542699194e+02,  2.0153708403353160716e-15},
  { 2.9132395009427028754e+02,  2.0029887665268810448e-14},
  { 2.9576660135076065217e+02, -2.8145997752905985118e-14},
  { 3.0022094864701415418e+02, -2.2421874778322229279e-14},
  { 3.0468685676566872189e+02, -6.4188626335518129545e-15},
  { 3.0916419358014690033e+02,  2.1613992872737720743e-14},
  { 3.1365282994987904885e+02,  1.2928984656414700999e-14},
  { 3.1815263962020929966e+02,  2.7185427014428510713e-14},
  { 3.2266349912672615119e+02,  2.5701792544522621036e-14},
  { 3.2718528770377520232e+02,  1.4883746199340294173e-14},
  { 3.3171788719692847280e+02,  3.3427014316894747085e-16},
  { 3.3626118197919845443e+02,  2.2606240796263950960e-14},
  { 3.4081505887079902095e+02, -3.0824154418584960781e-15},
  { 3.4537940706226686416e+02, -1.0051581128654682651e-14},
  { 3.4995411804077025408e+02, -1.7153763378924380938e-14},
  { 3.5453908551944078908e+02,  1.9764280928464860895e-14},
  { 3.5913420536957539753e+02,  1.2498389587699469416e-15},
  { 3.6373937555556346979e+02,  2.0355612947496250654e-14},
  { 3.6835449607240474279e+02,  6.8047236300469703549e-15},
  { 3.7297946888568901613e+02,  4.5427704126937566452e-15},
  { 3.7761419787391866976e+02, -1.3309881730125543788e-14},
  { 3.8225858877306001204e+02,  1.7066220149046898877e-14},
  { 3.8691254912321755910e+02, -6.6172288712915584272e-15},
  { 3.9157598821732960914e+02,  1.0483437695430516334e-14},
  { 3.9624881705179154778e+02, -2.1984000934124987299e-14},
  { 4.0093094827891576415e+02, -1.8657045236097604286e-14},
  { 4.0562229616114490227e+02, -1.3078724530527678942e-14},
  { 4.1032277652693733216e+02, -2.6743834869510098898e-14},
  { 4.1503230672824963676e+02,  2.7938237970287591604e-15},
  { 4.1975080559954471937e+02,  1.4728185592688318861e-14},
  { 4.2447819341825709216e+02, -1.7494428961102440175e-14},
  { 4.2921439186665156740e+02,  2.7308470015966966160e-15},
  { 4.3395932399501481314e+02,  7.0582536719504341250e-15},
  { 4.3871291418612116786e+02,  1.6981110787446912618e-14},
  { 4.4347508812091894015e+02,  8.0991519063490743489e-16},
  { 4.4824577274538461324e+02, -7.5260347471342305159e-15},
  { 4.5302489623849612599e+02,  9.1155395179683028364e-15},
  { 4.5781238798127816381e+02,  1.7291826846708180241e-14},
  { 4.6260817852687489449e+02,  2.7699251040825051093e-14},
  { 4.6741219957160819831e+02, -1.9567194605746345946e-14},
  { 4.7222438392698057896e+02,  1.7277261394956819310e-14},
  { 4.7704466549258563646e+02, -3.3523696904884624138e-15},
  { 4.8187297922988796017e+02, -2.5944561547034108804e-14},
  { 4.8670926113683941594e+02, -3.7140783665509818851e-15},
  { 4.9155344822329800536e+02, -1.8606661713197668769e-15},
  { 4.9640547848721763557e+02, -1.4909986887389890098e-14},
  { 5.0126529089157929775e+02, -4.9698975572564981668e-15},
  { 5.0613282534203489149e+02, -1.6291761523356007577e-14},
  { 5.1100802266523601247e+02,  1.4272158439089545512e-14},
  { 5.1589082458782240792e+02, -1.0322438559879302724e-14},
  { 5.2078117371604412256e+02,  2.8806364885275583495e-14},
  { 5.2567901351599505233e+02,  1.0402121035255553969e-14},
  { 5.3058428829443346331e+02,  2.8874529748818994254e-14},
  { 5.3549694318016952366e+02,  2.0532421502191848234e-14},
  { 5.4041692410599762297e+02,  4.6131459616760924960e-14},
  { 5.4534417779115483427e+02,  3.9529379627128233400e-14},
  { 5.5027865172428550977e+02,  5.5781769411128246342e-14},
  { 5.5522029414689484383e+02,  2.6018165361515851207e-14},
  { 5.6016905403727298562e+02,  5.2506102370160487189e-14},
  { 5.6512488109487435395e+02, -5.5087108956392259896e-14},
  { 5.7008772572513419163e+02,  1.4506980471415534004e-14},
  { 5.7505753902471019501e+02,  1.1753781140168979974e-14},
  { 5.8003427276713080118e+02, -2.0021327589606492820e-14},
  { 5.8501787938883910556e+02,  1.2043285070600097737e-14},
  { 5.9000831197561785757e+02, -3.6661231272690423630e-15},
  { 5.9500552424938200602e+02, -3.7052038350812133456e-14},
  { 6.0000947055532742525e+02,  2.8627630708621120029e-15},
  { 6.0502010584942365767e+02,  2.6184503415138050230e-14},
  { 6.1003738568623862193e+02, -1.3742309432483183413e-14},
  { 6.1506126620708494102e+02, -5.6440607359865996168e-14},
  { 6.2009170412847731768e+02,  2.3568246391384379285e-15},
  { 6.2512865673089095253e+02, -3.3342448773746644282e-15},
  { 6.3017208184781020464e+02, -8.8184229472838861186e-15},
  { 6.3522193785505976393e+02, -3.1069773573369826175e-14},
  { 6.4027818366040798992e+02,  5.1001997173422213274e-14},
  { 6.4534077869343502698e+02, -1.9255814291773831866e-14},
  { 6.5040968289565523719e+02,  2.0608186304978750433e-15},
  { 6.5548485671088906201e+02,  4.1569509299638154056e-15},
  { 6.6056626107587351271e+02,  1.6458398422786870119e-14},
  { 6.6565385741110594608e+02, -3.2842098384444132746e-14},
  { 6.7074760761191271285e+02, -3.7278289669647126945e-14},
  { 6.7584747403973688051e+02, -6.5139886982178939986e-15},
  { 6.8095341951363741373e+02,  4.0880286179157333063e-14},
  { 6.8606540730199401423e+02, -1.6390723998835485797e-14},
  { 6.9118340111441079898e+02, -4.6032539552775565324e-14},
  { 6.9630736509381404176e+02, -2.9882211401566375448e-14},
  { 7.0143726380873704329e+02,  4.2055124399124326631e-14},
  { 7.0657306224578735510e+02, -7.9849496927747738301e-15},
  { 7.1171472580228999050e+02,  1.6456292625473995897e-14},
  { 7.1686222027910343968e+02,  2.0316226926333451574e-14},
  { 7.2201551187360121276e+02,  2.6187160268978199580e-14},
  { 7.2717456717281572764e+02,  4.0325983740024313376e-14},
  { 7.3233935314673931316e+02, -3.1135586154860034624e-14},
  { 7.3750983714177743877e+02, -4.9669193222948042432e-15},
  { 7.4268598687435121519e+02,  4.7756515956293703443e-14},
  { 7.4786777042464336773e+02, -1.9636099794561631333e-14},
  { 7.5305515623048415819e+02, -5.5093579963689320400e-14},
  { 7.5824811308137429933e+02,  1.4142092581408200867e-14},
  { 7.6344661011264008721e+02,  5.2006783380256467323e-14},
  { 7.6865061679971688591e+02,  4.8658960236626277484e-14},
  { 7.7386010295255834990e+02,  5.6084809456045125924e-15},
  { 7.7907503871016729136e+02,  4.9769664787505897282e-14},
  { 7.8429539453524569126e+02, -2.5314216281309145887e-14},
  { 7.8952114120895885208e+02,  1.5115537849603218618e-14},
  { 7.9475224982581346467e+02, -1.0852001709196594192e-14},
  { 7.9998869178864345031e+02, -4.7287638590633906966e-14},
  { 8.0523043880370300940e+02,  3.6002229806206865392e-14},
  { 8.1047746287586357994e+02, -4.8394033966447898650e-14},
  { 8.1572973630391015831e+02,  3.1047649602296439976e-15},
  { 8.2098723167593789185e+02,  5.1114829460517171739e-14},
  { 8.2624992186484280410e+02,  2.4419303684157594321e-14},
  { 8.3151778002390619804e+02, -4.1390486100644422154e-14},
  { 8.3679077958246989510e+02,  8.3470351496137130874e-15},
  { 8.4206889424170037728e+02,  4.3400890783278654759e-14},
  { 8.4735209797043842173e+02, -1.2539129323500232650e-14},
  { 8.5264036500113297734e+02, -3.2917759966919451893e-14},
  { 8.5793366982585746428e+02, -2.7456778414557014287e-14},
  { 8.6323198719240542687e+02,  4.6628225761253004528e-14},
  { 8.6853529210046451681e+02,  3.2439948266782605572e-14},
  { 8.7384355979786573698e+02,  1.7024870275609020377e-14},
  { 8.7915676577690749127e+02,  5.0067221229110876502e-14},
  { 8.8447488577075171179e+02,  4.5941079338033365786e-14},
  { 8.8979789574989013090e+02,  3.5007414464190639557e-14},
  { 8.9512577191867978854e+02, -4.1552845213160979064e-14},
  { 9.0045849071194516000e+02, -4.3938958847244038476e-14},
  { 9.0579602879164644946e+02, -1.5427045641795900705e-14},
  { 9.1113836304361120710e+02,  3.7936062840264655799e-14},
  { 9.1648547057432870133e+02,  1.2388425146201894389e-14},
  { 9.2183732870780477242e+02,  7.7953594451063185545e-15},
  { 9.2719391498247682648e+02, -3.3811581522703927238e-14},
  { 9.3255520714818624128e+02, -2.3501988772091012436e-14},
  { 9.3792118316320807025e+02, -9.9019721456280313547e-16},
  { 9.4329182119133577089e+02, -3.8828228758975301703e-14},
  { 9.4866709959901993443e+02, -3.7365461734672225163e-14},
  { 9.5404699695256033465e+02,  2.1967819878695887981e-14},
  { 9.5943149201534947679e+02, -3.1161466950266501707e-14},
  { 9.6482056374516594133e+02,  5.1191120654897209981e-15},
  { 9.7021419129151831839e+02, -1.0410186947217000697e-14},
  { 9.7561235399303609483e+02, -3.4034480133133905328e-14},
  { 9.8101503137490828976e+02,  5.0482408720488663100e-14},
  { 9.8642220314636847434e+02, -1.5331743595366225784e-14},
  { 9.9183384919822344727e+02,  5.1581848387377046652e-14},
  { 9.9724994960042795356e+02, -3.4567282634413414685e-14},
  { 1.0026704845997002167e+03, -1.1810667001120463745e-14},
  { 1.0080954346171815814e+03,  2.6125277130975405712e-14},
  { 1.0135247802461360607e+03, -1.2357054440857773318e-14},
  { 1.0189585022496902411e+03,  4.6845014892965792221e-14},
  { 1.0243965815586134340e+03,  4.9308232889595227201e-14},
  { 1.0298389992691352290e+03,  4.7833129728631755020e-14},
  { 1.0352857366408015878e+03, -9.5162069040408270179e-16},
  { 1.0407367750943672036e+03,  8.3831400150598552722e-14},
  { 1.0461920962097249230e+03,  6.5796015880062249546e-14},
  { 1.0516516817238691601e+03, -1.2337586846647020826e-14},
  { 1.0571155135288947804e+03, -2.2537725458140462686e-14},
  { 1.0625835736700298639e+03,  2.5107267153829157218e-14},
  { 1.0680558443437014375e+03, -7.3802184976408644065e-14},
  { 1.0735323078956328118e+03,  6.2624111559301713342e-14},
  { 1.0790129468189747968e+03,  6.8936210214317284275e-14},
  { 1.0844977437524655670e+03, -4.6311548062238322692e-14},
  { 1.0899866814786221312e+03,  7.5918525179770268271e-14},
  { 1.0954797429219627247e+03,  3.0833994728088903865e-14},
  { 1.1009769111472560326e+03, -7.5140058453014808206e-14},
  { 1.1064781693578006525e+03,  3.1902669056994204030e-14},
  { 1.1119835008937329803e+03,  6.6882936383973062813e-14},
  { 1.1174928892303610155e+03,  8.9024535173005024901e-15},
  { 1.1230063179765259065e+03,  1.0009547147974694482e-13},
  { 1.1285237708729907808e+03, -6.6646410064312271630e-14},
  { 1.1340452317908529949e+03, -3.4259371733774208834e-14},
  { 1.1395706847299848050e+03, -6.0499280823710872181e-14},
  { 1.1451001138174960943e+03,  7.3553851404663992204e-14},
  { 1.1506335033062237017e+03, -1.3597233787223318058e-14},
  { 1.1561708375732421246e+03,  1.0007882745589962179e-13},
  { 1.1617121011184005965e+03,  5.4254719719550151903e-14}
};

// raise to 15 if we ever want sin_taylor() or cos_taylor()
CONSTI32(_ddr_n_inv_fact, 6);

// "_offset3" since this starts at 1 / (3!), not 1 / (0!)
static const dd_real  _ddr_inv_fact_offset3[_ddr_n_inv_fact] = {
  { 1.66666666666666657e-01,  9.25185853854297066e-18},
  { 4.16666666666666644e-02,  2.31296463463574266e-18},
  { 8.33333333333333322e-03,  1.15648231731787138e-19},
  { 1.38888888888888894e-03, -5.30054395437357706e-20},
  { 1.98412698412698413e-04,  1.72095582934207053e-22},
  { 2.48015873015873016e-05,  2.15119478667758816e-23}
};

dd_real ddr_exp(const dd_real a) {
  // Strategy: We first reduce the size of x by noting that
  //
  //   exp(kr + m * log(2)) = 2^m * exp(r)^k
  //
  // where m and k are integers.  By choosing m appropriately
  // we can make |kr| <= log(2) / 2 = 0.347.  Then exp(r) is
  // evaluated using the familiar Taylor series.  Reducing the
  // argument substantially speeds up the convergence.

  const double k = 512.0;
  const double inv_k = 1.0 / k;

  if (a.x[0] <= -709.0) {
    return ddr_maked(0.0);
  }

  if (a.x[0] >= 709.0) {
    return ddr_make(INFINITY_D, INFINITY_D);
  }

  if (ddr_is_zero(a)) {
    return ddr_maked(1.0);
  }

  if (ddr_is_one(a)) {
    return _ddr_e;
  }

  double m = floor(a.x[0] / _ddr_log2.x[0] + 0.5);
  dd_real r = ddr_mul_pwr2(ddr_sub(a, ddr_muld(_ddr_log2, m)), inv_k);
  dd_real s, t, p;

  p = ddr_sqr(r);
  s = ddr_add(r, ddr_mul_pwr2(p, 0.5));
  p = ddr_mul(p, r);
  t = ddr_mul(p, _ddr_inv_fact_offset3[0]);
  int32_t i = 0;
  do {
    s = ddr_add(s, t);
    p = ddr_mul(p, r);
    ++i;
    t = ddr_mul(p, _ddr_inv_fact_offset3[i]);
  } while ((fabs(t.x[0]) > inv_k * _ddr_eps) && (i < 5));

  s = ddr_add(s, t);

  s = ddr_add(ddr_mul_pwr2(s, 2.0), ddr_sqr(s));
  s = ddr_add(ddr_mul_pwr2(s, 2.0), ddr_sqr(s));
  s = ddr_add(ddr_mul_pwr2(s, 2.0), ddr_sqr(s));
  s = ddr_add(ddr_mul_pwr2(s, 2.0), ddr_sqr(s));
  s = ddr_add(ddr_mul_pwr2(s, 2.0), ddr_sqr(s));
  s = ddr_add(ddr_mul_pwr2(s, 2.0), ddr_sqr(s));
  s = ddr_add(ddr_mul_pwr2(s, 2.0), ddr_sqr(s));
  s = ddr_add(ddr_mul_pwr2(s, 2.0), ddr_sqr(s));
  s = ddr_add(ddr_mul_pwr2(s, 2.0), ddr_sqr(s));
  s = ddr_addd(s, 1.0);

  return ddr_ldexp(s, S_CAST(int32_t, m));
}

dd_real ddr_log(const dd_real a) {
  // Strategy.  The Taylor series for log converges much more
  // slowly than that of exp, due to the lack of the factorial
  // term in the denominator.  Hence this routine instead tries
  // to determine the root of the function
  //
  //   f(x) = exp(x) - a
  //
  // using Newton iteration.  The iteration is given by
  //
  //   x' = x - f(x)/f'(x)
  //      = x - (1 - a * exp(-x))
  //      = x + a * exp(-x) - 1.
  //
  // Only one iteration is needed, since Newton's iteration
  // approximately doubles the number of digits per iteration.

  if (ddr_is_one(a)) {
    return ddr_maked(0.0);
  }

  if (a.x[0] <= 0.0) {
    return ddr_make(0.0 / 0.0, 0.0 / 0.0);
  }

  dd_real x;
  x.x[0] = log(a.x[0]);  // Initial approximation
  x.x[1] = 0.0;

  x = ddr_subd(ddr_add(x, ddr_mul(a, ddr_exp(ddr_negate(x)))), 1.0);
  return x;
}

dd_real ddr_lfact(double xx) {
  // It can be shown from the Euler-Maclaurin approximation that
  //   ln (n!) = n ln n - n + 0.5 ln (2*pi*n) + n^{-1}/12 - n^{-3}/360
  //             + n^{-5}/1260 - n^{-7}/1680 + n^{-9}/1188 - 691*n^{-11}/360360
  //             + ...
  // and that if we truncate the calculation, the error is bounded by the
  // magnitude of the first omitted term for non-tiny n.
  //
  // The terms decrease more quickly for larger n, so it is reasonable to use a
  // lookup table for small n.  We stop the lookup table at 256 since the
  // -691*n^{-11}/360360 term (which we omit) becomes too small to be relevant
  // to double-double accuracy around that point.
  if (xx < 256) {
    return _ddr_ln_fact[S_CAST(int32_t, xx)];
  }
  // Since we're summing a fixed number of terms, start with the smaller terms
  // to reduce rounding error.
  const dd_real invn = ddr_divd(ddr_maked(1.0), xx);
  const dd_real invn2 = ddr_sqr(invn);
  dd_real sum =
    ddr_mul(invn,
            ddr_add(_ddr_12th,
                    ddr_mul(invn2,
                            ddr_add(_ddr_neg_360th,
                                    ddr_mul(invn2,
                                            ddr_add(_ddr_1260th,
                                                    ddr_mul(invn2,
                                                            ddr_add(_ddr_neg_1680th,
                                                                    ddr_mul(invn2, _ddr_1188th)))))))));
  //   n ln n - n + 0.5 ln (2*pi*n)
  // = n ln n - n + 0.5 ln 2*pi + 0.5 ln n
  // = (n + 0.5) ln n - n + 0.5 ln 2*pi
  sum = ddr_add(sum, _ddr_half_log_2pi);
  sum = ddr_subd(sum, xx);
  const dd_real logn = ddr_log(ddr_maked(xx));
  return ddr_add(sum, ddr_muld(logn, xx + 0.5));
}

// Bignum factorial helper functions and constants.
// See GMP mpz/bin_uiui.c.
#if defined(__LP64__) && !defined(_WIN32)
// assume first term < 2^32 for now and make this comparison compile to no-op
static const mp_limb_t kFallingFactorialInitMul2Max = 0xffffffffffffffffLLU;

static const mp_limb_t kFallingFactorialInitMul3Max = 0x32cbfe;
static const mp_limb_t kFallingFactorialInitMul4Max = 0x16a0b;
static const mp_limb_t kFallingFactorialInitMul5Max = 0x24c4;
static const mp_limb_t kFallingFactorialInitMul6Max = 0xa16;
static const mp_limb_t kFallingFactorialInitMul7Max = 0x34b;
static const mp_limb_t kFallingFactorialInitMul8Max = 0x1b2;
#else
static const mp_limb_t kFallingFactorialInitMul2Max = 0x16a0a;
static const mp_limb_t kFallingFactorialInitMul3Max = 0x801;
static const mp_limb_t kFallingFactorialInitMul4Max = 0x16b;
static const mp_limb_t kFallingFactorialInitMul5Max = 0x73;
static const mp_limb_t kFallingFactorialInitMul6Max = 0x42;
static const mp_limb_t kFallingFactorialInitMul7Max = 0x26;
// 0x1e doesn't work when ct is just 7
static const mp_limb_t kFallingFactorialInitMul8Max = 0x1a;
#endif

// mul2..mul8 remove a constant (1, 1, 2, 2, 4, 4, 6) number of factors of 2
// along the way; those can be restored with a single left-shift at the end.
static inline mp_limb_t mul2(mp_limb_t m) {
  /* We need to shift before multiplying, to avoid an overflow. */
  mp_limb_t m01 = (m | 1) * ((m + 1) >> 1);
  return m01;
}

static inline mp_limb_t mul3(mp_limb_t m) {
  mp_limb_t m01 = (m + 0) * (m + 1) >> 1;
  mp_limb_t m2 = (m + 2);
  return m01 * m2;
}

static inline mp_limb_t mul4(mp_limb_t m) {
  mp_limb_t m03 = (m + 0) * (m + 3) >> 1;
  return m03 * (m03 + 1); /* mul2 (m03) ? */
}

static inline mp_limb_t mul5(mp_limb_t m) {
  mp_limb_t m03 = (m + 0) * (m + 3) >> 1;
  mp_limb_t m034 = m03 * (m + 4);
  return (m03 + 1) * m034;
}

static inline mp_limb_t mul6(mp_limb_t m) {
  mp_limb_t m05 = (m + 0) * (m + 5);
  mp_limb_t m1234 = (m05 + 5) * (m05 + 5) >> 3;
  return m1234 * (m05 >> 1);
}

static mp_limb_t mul7(mp_limb_t m) {
  mp_limb_t m05 = (m + 0) * (m + 5);
  mp_limb_t m1234 = (m05 + 5) * (m05 + 5) >> 3;
  mp_limb_t m056 = m05 * (m + 6) >> 1;
  return m1234 * m056;
}

static mp_limb_t mul8(mp_limb_t m) {
  mp_limb_t m07 = (m + 0) * (m + 7);
  mp_limb_t m0257 = m07 * (m07 + 10) >> 3;
  mp_limb_t m1346 = m07 + 9 + m0257;
  return m0257 * m1346;
}

// Precondition: main[0], main[1], main[2], ..., main[num_ct - 1] are the
//               positive numbers to multiply.  main and wkspace are both
//               allocated up to at least [num_ct rounded up to power of 2]
//               mp_limb_ts.  No interior number is larger than max(main[0],
//               main[num_ct - 1]).
//
// Postcondition: main has the multiplied result.  Return value is final limb
//                count.
uint32_t bottom_up_multiply(uint32_t num_ct, mp_limb_t* main, mp_limb_t* wkspace) {
  mp_limb_t* src = main;
  mp_limb_t* dst = wkspace;
  uint32_t input_limb_ct = 1;
  uint32_t input_stride = 1;
  while (num_ct > 1) {
    mp_limb_t* src_iter = src;
    mp_limb_t* dst_iter = dst;
    uint32_t next_num_ct = num_ct / 2;
    uint32_t next_stride = input_limb_ct * 2;
    for (uint32_t num_idx = 0; num_idx != next_num_ct; ++num_idx) {
      mp_limb_t* src2 = &(src_iter[input_stride]);
      mpn_mul_n(dst_iter, src_iter, src2, input_limb_ct);
      src_iter = &(src2[input_stride]);
      dst_iter = &(dst_iter[next_stride]);
    }
    if (num_ct % 2) {
      memcpy(dst_iter, src_iter, input_limb_ct * sizeof(mp_limb_t));
      memset(&(dst_iter[input_limb_ct]), 0, (next_stride - input_limb_ct) * sizeof(mp_limb_t));
      dst_iter = &(dst_iter[next_stride]);
      ++next_num_ct;
    }
    mp_limb_t* swap_ptr = src;
    src = dst;
    dst = swap_ptr;
    input_limb_ct = next_stride;
    // If both the first and the last numbers have top limb 0, all other
    // numbers must also have top limb 0 due to the last part of the
    // precondition.
    dst_iter = &(dst_iter[-S_CAST(intptr_t, next_stride)]);
    while ((src[input_limb_ct - 1] == 0) && (dst_iter[input_limb_ct - 1] == 0)) {
      --input_limb_ct;
    }
    input_stride = next_stride;
    num_ct = next_num_ct;
  }
  if (src != main) {
    memcpy(main, src, input_limb_ct * sizeof(mp_limb_t));
  }
  return input_limb_ct;
}

// Assumes result and wkspace are each large enough to fit [ct rounded up to
// power of 2] 32-bit integers.
//
// Return value is size of final left-shift needed to produce correct
// falling-factorial.
int32_t falling_factorial(mp_limb_t top, uint32_t ct, mp_limb_t* result, uint32_t* result_limb_ctp, mp_limb_t* wkspace) {
  const mp_limb_t stop = top - ct;
  uint32_t final_lshift = 0;
  uint32_t num_ct;
  if (top > kFallingFactorialInitMul5Max) {
    if (top > kFallingFactorialInitMul3Max) {
      if (top > kFallingFactorialInitMul2Max) {
        num_ct = ct;
        for (uint32_t num_idx = 0; num_idx != num_ct; ++num_idx) {
          result[num_idx] = top;
          --top;
        }
      } else {
        num_ct = ct / 2;
        top -= 1;
        for (uint32_t num_idx = 0; num_idx != num_ct; ++num_idx) {
          result[num_idx] = mul2(top);
          top -= 2;
        }
        final_lshift = 1 * num_ct;
        top += 1;
      }
    } else {
      if (top > kFallingFactorialInitMul4Max) {
        num_ct = ct / 3;
        top -= 2;
        for (uint32_t num_idx = 0; num_idx != num_ct; ++num_idx) {
          result[num_idx] = mul3(top);
          top -= 3;
        }
        final_lshift = 1 * num_ct;
        top += 2;
      } else {
        num_ct = ct / 4;
        top -= 3;
        for (uint32_t num_idx = 0; num_idx != num_ct; ++num_idx) {
          result[num_idx] = mul4(top);
          top -= 4;
        }
        final_lshift = 2 * num_ct;
        top += 3;
      }
    }
  } else {
    if (top > kFallingFactorialInitMul7Max) {
      if (top > kFallingFactorialInitMul6Max) {
        num_ct = ct / 5;
        top -= 4;
        for (uint32_t num_idx = 0; num_idx != num_ct; ++num_idx) {
          result[num_idx] = mul5(top);
          top -= 5;
        }
        final_lshift = 2 * num_ct;
        top += 4;
      } else {
        num_ct = ct / 6;
        top -= 5;
        for (uint32_t num_idx = 0; num_idx != num_ct; ++num_idx) {
          result[num_idx] = mul6(top);
          top -= 6;
        }
        final_lshift = 4 * num_ct;
        top += 5;
      }
    } else {
      if (top > kFallingFactorialInitMul8Max) {
        num_ct = ct / 7;
        top -= 6;
        for (uint32_t num_idx = 0; num_idx != num_ct; ++num_idx) {
          result[num_idx] = mul7(top);
          top -= 7;
        }
        final_lshift = 4 * num_ct;
        top += 6;
      } else {
        num_ct = ct / 8;
        top -= 7;
        for (uint32_t num_idx = 0; num_idx != num_ct; ++num_idx) {
          result[num_idx] = mul8(top);
          top -= 8;
        }
        final_lshift = 6 * num_ct;
        top += 7;
      }
    }
  }
  if (top > stop) {
    mp_limb_t last_term = top;
    while (--top > stop) {
      last_term *= top;
    }
    result[num_ct++] = last_term;
  }
  *result_limb_ctp = bottom_up_multiply(num_ct, result, wkspace);
  return final_lshift;
}

void lshift_multilimb(uint64_t lshift_ct, mp_limb_t* num, uint32_t* num_limb_ctp) {
  const uintptr_t limb_shift_ct = lshift_ct / mp_bits_per_limb;
  const uint32_t lshift_rem = lshift_ct % mp_bits_per_limb;
  mp_limb_t* num_write = &(num[limb_shift_ct]);
  uint32_t num_limb_ct = *num_limb_ctp;
  if (lshift_rem) {
    const mp_limb_t new_top_limb = mpn_lshift(num_write, num, num_limb_ct, lshift_rem);
    if (new_top_limb) {
      num_write[num_limb_ct++] = new_top_limb;
    }
  } else {
    memmove(num_write, num, num_limb_ct * sizeof(mp_limb_t));
  }
#ifdef __LP64__
  assert(num_limb_ct + limb_shift_ct <= 0xffffffffU);
#endif
  memset(num, 0, limb_shift_ct * sizeof(mp_limb_t));
  *num_limb_ctp = num_limb_ct + limb_shift_ct;
}

// Preconditions:
// - numer_factorial_args[] and denom_factorial_args[] are
//   not-necessarily-sorted lists of length ffac_ct, describing a quotient of
//   factorial-products.  (If one list is longer than the other, just pad the
//   other with zeroes.)
// - pow2 is a power of 2 to multiply the quotient by at the end.
// - starting_lnprobv_ddr is either initialized to
//     log(2^numer_pow2 / (numer_factorial_args[0]! ... numer_factorial_args[ffac_ct-1]!))
//   or it has x[0] initialized to DBL_MAX to indicate that the calculation
//   hasn't happened.  In the latter case, it may be set to the former value
//   if that is needed.
//
// This function errors out iff memory allocation fails.
//
// Postconditions on success:
// - *cmp_resultp is set to a positive value if the fraction > 1, a negative
//   value if the fraction < 1, and zero if it's exactly 1.
// - *dbl_ptr is the double representation of the fraction, error limited to
//   1-2 ulps.
// - numer_factorial_args[] and denom_factorial_args[] are sorted in
//   nondecreasing order.
//
// This could take a precomputed ddr_lfact table as an additional pair of
// parameters, but I don't think that makes much of a difference.
BoolErr CompareFactorialProducts(uint32_t ffac_ct, int64_t pow2, int64_t numer_pow2, uint32_t* numer_factorial_args, uint32_t* denom_factorial_args, dd_real* starting_lnprobv_ddr_ptr, mp_limb_t** gmp_wkspacep, uintptr_t* gmp_wkspace_limb_ctp, intptr_t* cmp_resultp, double* dbl_ptr) {
  // 1. Sort numer_factorial_args and denom_factorial_args.  (This has the
  //    effect of cancelling out matching terms.)
  // 2. Iterate through numer_factorial_args[] and denom_factorial_args[] to
  //    determine bignum calculation size.
  // 3. If bignum calculation size is large enough, perform the comparison with
  //    dd_reals first, returning a result unless the log-likelihoods are
  //    within ffac_ct * 2^{-60} of each other.
  // 4. Reallocate workspace if necessary, and iterate properly through
  //    numer_factorial_args[] and denom_factorial_args[].  When
  //    numer_factorial_args[k] > denom_factorial_args[k], multiply the
  //    numerator by the corresponding falling-factorial; and vice versa.
  // 5. Perform final left-shift of numerator or denominator (necessary since
  //    we allow a 2^pow2 term in the expression, and we pull out some factors
  //    of 2 while accumulating the falling-factorial products).
  // 6. Compare.
  //
  // Possible improvement: find appropriate spots to use mpn_sqr().
  STD_SORT(ffac_ct, u32cmp, numer_factorial_args);
  STD_SORT(ffac_ct, u32cmp, denom_factorial_args);

  uintptr_t numer_term_ct = 0;
  uintptr_t denom_term_ct = 0;
  uint32_t max_ffac_size = 0;
  for (uint32_t ffac_idx = 0; ffac_idx != ffac_ct; ++ffac_idx) {
    const uint32_t numer_factorial_arg = numer_factorial_args[ffac_idx];
    const uint32_t denom_factorial_arg = denom_factorial_args[ffac_idx];
    uint32_t ffac_size;
    if (numer_factorial_arg > denom_factorial_arg) {
      ffac_size = numer_factorial_arg - denom_factorial_arg;
      numer_term_ct += ffac_size;
    } else {
      ffac_size = denom_factorial_arg - numer_factorial_arg;
      denom_term_ct += ffac_size;
    }
    if (max_ffac_size < ffac_size) {
      max_ffac_size = ffac_size;
    }
  }
  if (max_ffac_size == 0) {
    // All factorials cancel out.
    *cmp_resultp = pow2;
    *dbl_ptr = kLn2 * S_CAST(double, pow2);
    return 0;
  }
  if (numer_term_ct + denom_term_ct > 256) {
    dd_real starting_lnprobv_ddr = *starting_lnprobv_ddr_ptr;
    if (starting_lnprobv_ddr.x[0] == DBL_MAX) {
      starting_lnprobv_ddr = ddr_muld(_ddr_log2, S_CAST(double, numer_pow2));
      for (uint32_t ffac_idx = 0; ffac_idx != ffac_ct; ++ffac_idx) {
        starting_lnprobv_ddr = ddr_sub(starting_lnprobv_ddr, ddr_lfact(u31tod(numer_factorial_args[ffac_idx])));
      }
      *starting_lnprobv_ddr_ptr = starting_lnprobv_ddr;
    }
    dd_real lnprobv_ddr = ddr_muld(_ddr_log2, S_CAST(double, numer_pow2 + pow2));
    for (uint32_t ffac_idx = 0; ffac_idx != ffac_ct; ++ffac_idx) {
      lnprobv_ddr = ddr_sub(lnprobv_ddr, ddr_lfact(u31tod(denom_factorial_args[ffac_idx])));
    }
    const double lnprob_diff = ddr_sub(lnprobv_ddr, starting_lnprobv_ddr).x[0];
    // ddr_lfact() result has >= 106 bits of precision, should be accurate to
    // 96+ bits (mostly dependent on log).
    // log((2^31)!) is less than 2^36, so we should have at least 60 bits past
    // the decimal point.
    const double epsilon = k2m60 * u31tod(ffac_ct);
    if (lnprob_diff > epsilon) {
      *cmp_resultp = 1;
      *dbl_ptr = exp(lnprob_diff);
      return 0;
    }
    if (lnprob_diff < -epsilon) {
      *cmp_resultp = -1;
      *dbl_ptr = exp(lnprob_diff);
      return 0;
    }
  }

  // numer: Usually CeilPow2(numer_term_ct) limbs in 32-bit limb case.  Also
  //        ensure this is at least numer_term_ct + max(1, ceil(pow2 / 32)), to
  //        cover mpn_mul() and lshift_multilimb() edge cases.
  //        Usually DivUp(CeilPow2(numer_term_ct), 2) limbs in 64-bit limb
  //        case.  Also ensure this is at least DivUp(numer_term_ct, 2) +
  //        max(1, ceil(pow2 / 64)) to cover edge cases.
  // denom: Similar to above.
  // new_ffac: CeilPow2(max_ffac_size) in 32-bit limb case, DivUp(., 2) in
  //           64-bit case.
  // generic_wkspace: Used as intermediate buffer for falling-factorial
  //                  calculation, and temporary buffer to copy previous
  //                  numerator or denominator into before multiplication by
  //                  new_falling_factorial.  Safe to make this the larger of
  //                  the numer and denom sizes; could tighten this later.
  const uintptr_t numer_bound1 = numer_term_ct? DivUp(CeilPow2(numer_term_ct), kInt32PerLimb) : 0;
  uintptr_t numer_bound2 = 1;
  const uintptr_t denom_bound1 = denom_term_ct? DivUp(CeilPow2(denom_term_ct), kInt32PerLimb) : 0;
  uintptr_t denom_bound2 = 1;
  if (pow2 > 0) {
    numer_bound2 = DivUp(pow2, mp_bits_per_limb);
  } else if (pow2 < 0) {
    denom_bound2 = DivUp(-pow2, mp_bits_per_limb);
  }
  numer_bound2 += numer_term_ct;
  denom_bound2 += denom_term_ct;
  const uintptr_t numer_limb_req = MAXV(numer_bound1, numer_bound2);
  const uintptr_t denom_limb_req = MAXV(denom_bound1, denom_bound2);
  const uintptr_t new_ffac_limb_req = DivUp(CeilPow2(max_ffac_size), kInt32PerLimb);
  const uintptr_t generic_wkspace_limb_req = MAXV(numer_limb_req, denom_limb_req);
  const uintptr_t total_limb_req = numer_limb_req + denom_limb_req + new_ffac_limb_req + generic_wkspace_limb_req;
  if (total_limb_req > *gmp_wkspace_limb_ctp) {
    free_cond(*gmp_wkspacep);
    uintptr_t new_limb_ct = 2 * (*gmp_wkspace_limb_ctp);
    if (new_limb_ct < total_limb_req) {
      new_limb_ct = total_limb_req;
    }
    *gmp_wkspacep = S_CAST(mp_limb_t*, malloc(total_limb_req * sizeof(mp_limb_t)));
    if (unlikely(!(*gmp_wkspacep))) {
      return 1;
    }
    *gmp_wkspace_limb_ctp = new_limb_ct;
  }
  mp_limb_t* numer = *gmp_wkspacep;
  mp_limb_t* denom = &(numer[numer_limb_req]);
  mp_limb_t* new_ffac = &(denom[denom_limb_req]);
  mp_limb_t* generic_wkspace = &(new_ffac[new_ffac_limb_req]);

  uint32_t numer_limb_ct = 1;
  uint32_t denom_limb_ct = 1;
  numer[0] = 1;
  denom[0] = 1;
  for (uint32_t ffac_idx = 0; ffac_idx != ffac_ct; ++ffac_idx) {
    const uint32_t numer_factorial_arg = numer_factorial_args[ffac_idx];
    const uint32_t denom_factorial_arg = denom_factorial_args[ffac_idx];
    if (numer_factorial_arg == denom_factorial_arg) {
      continue;
    }
    mp_limb_t* main;
    uint32_t* main_limb_ct_ptr;
    uint32_t top;
    uint32_t ct;
    int32_t sign;
    if (numer_factorial_arg > denom_factorial_arg) {
      main = numer;
      main_limb_ct_ptr = &numer_limb_ct;
      top = numer_factorial_arg;
      ct = top - denom_factorial_arg;
      sign = 1;
    } else {
      main = denom;
      main_limb_ct_ptr = &denom_limb_ct;
      top = denom_factorial_arg;
      ct = top - numer_factorial_arg;
      sign = -1;
    }
    uint32_t new_ffac_limb_ct;
    pow2 += sign * falling_factorial(top, ct, new_ffac, &new_ffac_limb_ct, generic_wkspace);
    const uint32_t main_limb_ct = *main_limb_ct_ptr;
    memcpy(generic_wkspace, main, main_limb_ct * sizeof(mp_limb_t));
    mp_limb_t msl;
    if (new_ffac_limb_ct >= main_limb_ct) {
      msl = mpn_mul(main, new_ffac, new_ffac_limb_ct, generic_wkspace, main_limb_ct);
    } else {
      msl = mpn_mul(main, generic_wkspace, main_limb_ct, new_ffac, new_ffac_limb_ct);
    }
    *main_limb_ct_ptr = main_limb_ct + new_ffac_limb_ct - (msl == 0);
  }
  if (pow2 > 0) {
    lshift_multilimb(pow2, numer, &numer_limb_ct);
  } else if (pow2 < 0) {
    lshift_multilimb(-pow2, denom, &denom_limb_ct);
  }
  double numer_d = numer[numer_limb_ct - 1];
  double denom_d = denom[denom_limb_ct - 1];
#ifdef __LP64__
  numer_d *= k2p64;
  if (numer_limb_ct > 1) {
    numer_d += numer[numer_limb_ct - 2];
  }
  denom_d *= k2p64;
  if (denom_limb_ct > 1) {
    denom_d += denom[denom_limb_ct - 2];
  }
#else
  numer_d *= 4294967296.0;
  if (numer_limb_ct > 1) {
    numer_d += numer[numer_limb_ct - 2];
    numer_d *= 4294967296.0;
    if (numer_limb_ct > 2) {
      numer_d += numer[numer_limb_ct - 3];
    }
  }
  denom_d *= 4294967296.0;
  if (denom_limb_ct > 1) {
    denom_d += denom[denom_limb_ct - 2];
    denom_d *= 4294967296.0;
    if (denom_limb_ct > 2) {
      denom_d += denom[denom_limb_ct - 3];
    }
  }
#endif
  double ratio = numer_d / denom_d;
  if (numer_limb_ct == denom_limb_ct) {
    *cmp_resultp = mpn_cmp(numer, denom, numer_limb_ct);
  } else {
    const int32_t limb_diff_ct = numer_limb_ct - denom_limb_ct;
    *cmp_resultp = S_CAST(int32_t, limb_diff_ct);
    ratio = scalbn(ratio, limb_diff_ct * mp_bits_per_limb);
  }
  *dbl_ptr = ratio;
  return 0;
}

#ifdef __cplusplus
}
#endif
