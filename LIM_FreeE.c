/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * LIM_FreeE.c
 *
 * Code generation for function 'LIM_FreeE'
 *
 */

/* Include files */
#include <math.h>
#include "rt_nonfinite.h"
#include <string.h>
#include "LIM_FreeE.h"
#include "LIM_FreeE_emxutil.h"

/* Function Declarations */
static double MutualInd(double ra, double rb, double h);
static void Mutual_dx_p2c_vec(double Coil_NumOfTurn, const emxArray_real_T
  *Coil_CoilMap, double Proj_NumOfElem, const emxArray_real_T *Proj_ElemMap,
  double Proj_Zp, emxArray_real_T *Mutual_dx);
static void RectCoil_RectCoil(double rp, double ac, double Num_lay, double lpn,
  double *obj_NumOfTurn, double *obj_SelfInd, double *obj_Resist,
  emxArray_real_T *obj_CoilMap, double *obj_Zc);
static void __anon_fcn(double CP_Lsw, double CP_LD, double Coil_NumOfTurn, const
  emxArray_real_T *Coil_CoilMap, double Proj_NumOfElem, const emxArray_real_T
  *Proj_ElemMap, double Proj_Mass, const emxArray_real_T *AA, const
  emxArray_real_T *bb, const emxArray_real_T *y, emxArray_real_T *varargout_1);
static void callODEFunctionNSM(double c_odeFcn_tunableEnvironment_f1_, double
  odeFcn_tunableEnvironment_f1_LD, double c_odeFcn_tunableEnvironment_f2_, const
  emxArray_real_T *d_odeFcn_tunableEnvironment_f2_, double
  c_odeFcn_tunableEnvironment_f3_, const emxArray_real_T
  *d_odeFcn_tunableEnvironment_f3_, double e_odeFcn_tunableEnvironment_f3_,
  const emxArray_real_T *odeFcn_tunableEnvironment_f4, const emxArray_real_T
  *odeFcn_tunableEnvironment_f5, const emxArray_real_T *y, emxArray_real_T *yp);
static double eps(double x);
static void invNxN(const emxArray_real_T *x, emxArray_real_T *y);
static void maxAbsThresh(const emxArray_real_T *a, const emxArray_real_T *b,
  emxArray_real_T *y);
static double norm(const emxArray_real_T *x);
static void ntrp45(const double t[3], double t0, const emxArray_real_T *b_y0,
                   double h, const emxArray_real_T *f, emxArray_real_T *y);
static void ode45(double ode_tunableEnvironment_f1_Lsw, double
                  ode_tunableEnvironment_f1_LD, double
                  c_ode_tunableEnvironment_f2_Num, const emxArray_real_T
                  *c_ode_tunableEnvironment_f2_Coi, double
                  c_ode_tunableEnvironment_f3_Num, const emxArray_real_T
                  *c_ode_tunableEnvironment_f3_Ele, double
                  ode_tunableEnvironment_f3_Mass, const emxArray_real_T
                  *ode_tunableEnvironment_f4, const emxArray_real_T
                  *ode_tunableEnvironment_f5, const double tspan[2], const
                  emxArray_real_T *b_y0, emxArray_real_T *varargout_1,
                  emxArray_real_T *varargout_2);
static void rdivide(const emxArray_real_T *x, const emxArray_real_T *y,
                    emxArray_real_T *z);
static double rt_powd_snf(double u0, double u1);

/* Function Definitions */
static double MutualInd(double ra, double rb, double h)
{
  double b0;
  double k;
  double s0;
  double a0;
  double i1;
  double w1;
  double a1;
  double x;

  /* Mutual inductance */
  b0 = ra + rb;
  k = sqrt(4.0 * ra * rb / (h * h + b0 * b0));
  s0 = k * k;
  if (s0 == 1.0) {
    a0 = rtInf;
    b0 = 1.0;
  } else if ((s0 >= 0.0) && (s0 <= 1.0)) {
    a0 = 1.0;
    b0 = sqrt(1.0 - s0);
    i1 = 0.0;
    w1 = 1.0;
    a1 = 1.0;
    while (w1 > 2.2204460492503131E-16) {
      a1 = (a0 + b0) / 2.0;
      x = a0 * b0;
      b0 = (a0 - b0) / 2.0;
      i1++;
      w1 = rt_powd_snf(2.0, i1) * b0 * b0;
      s0 += w1;
      a0 = a1;
      b0 = sqrt(x);
    }

    a0 = 3.1415926535897931 / (2.0 * a1);
    b0 = a0 * (1.0 - s0 / 2.0);
  } else {
    a0 = rtNaN;
    b0 = rtNaN;
  }

  return 1.2566370614359173E-6 * sqrt(ra * rb) * ((2.0 / k - k) * a0 - 2.0 / k *
    b0);
}

static void Mutual_dx_p2c_vec(double Coil_NumOfTurn, const emxArray_real_T
  *Coil_CoilMap, double Proj_NumOfElem, const emxArray_real_T *Proj_ElemMap,
  double Proj_Zp, emxArray_real_T *Mutual_dx)
{
  int j;
  int loop_ub;
  double h;
  double sum;
  double w1;
  int b_j;
  double theta;
  static const double We[20] = { -0.97390652851717174, -0.86506336668898454,
    -0.67940956829902444, -0.43339539412924721, -0.14887433898163122,
    0.14887433898163122, 0.43339539412924721, 0.67940956829902444,
    0.86506336668898454, 0.97390652851717174, 0.066671344308688138,
    0.14945134915058059, 0.21908636251598204, 0.26926671930999635,
    0.29552422471475287, 0.29552422471475287, 0.26926671930999635,
    0.21908636251598204, 0.14945134915058059, 0.066671344308688138 };

  /* dM/dxを出力すると、後の電磁力を算出するフェーズで結果を再利用できそうなので、検討しておいて。 */
  /* てか、そのようにしましょう。 */
  j = Mutual_dx->size[0];
  Mutual_dx->size[0] = (int)Proj_NumOfElem;
  emxEnsureCapacity_real_T1(Mutual_dx, j);
  loop_ub = (int)Proj_NumOfElem;
  for (j = 0; j < loop_ub; j++) {
    Mutual_dx->data[j] = 0.0;
  }

  for (j = 0; j < (int)Proj_NumOfElem; j++) {
    for (loop_ub = 0; loop_ub < (int)Coil_NumOfTurn; loop_ub++) {
      h = (Coil_CoilMap->data[loop_ub + Coil_CoilMap->size[0]] -
           Proj_ElemMap->data[j + Proj_ElemMap->size[0]]) - Proj_Zp;

      /* Gauuss lugandle */
      sum = 0.0;
      w1 = 2.0 * Coil_CoilMap->data[loop_ub] * Proj_ElemMap->data[j] / ((h * h +
        Coil_CoilMap->data[loop_ub] * Coil_CoilMap->data[loop_ub]) +
        Proj_ElemMap->data[j] * Proj_ElemMap->data[j]);
      for (b_j = 0; b_j < 10; b_j++) {
        theta = 1.5707963267948966 * We[b_j] + 1.5707963267948966;
        sum += We[10 + b_j] * (cos(theta) * ((1.0 - w1) * (1.0 - w1))) /
          rt_powd_snf(1.0 - w1 * cos(theta), 1.5);
      }

      /* Mdx=sum*(mu/(2^(3/2)))*(h/sqrt(a*b*w1))*(w1/(1-w1))^2*pi; */
      theta = w1 / (1.0 - w1);

      /* %%%積分範囲で0.5を掛けてみた　 */
      Mutual_dx->data[j] += sum * 4.4428829381583664E-7 * (h / sqrt
        (Coil_CoilMap->data[loop_ub] * Proj_ElemMap->data[j] * w1)) * (theta *
        theta) * 3.1415926535897931 * 0.5;
    }
  }
}

static void RectCoil_RectCoil(double rp, double ac, double Num_lay, double lpn,
  double *obj_NumOfTurn, double *obj_SelfInd, double *obj_Resist,
  emxArray_real_T *obj_CoilMap, double *obj_Zc)
{
  double rc;
  int i;
  int loop_ub;
  double Total;

  /*        function CoilMap = CoilGene(rp,ac,Num_lay,lpn) */
  /*           CoilMap = CoilGene(rp,ac,Num_lay,lpn); */
  /*        end */
  /* Coil data gene */
  /* 矩形断面コイルデータ作成 */
  /* zc_w=0;             %最下層コイルの最手前座標 */
  /* rp=20*1e-3;         %プロジェクタイル外径 */
  /* ac=1*1e-3;          %最下層コイル巻き線太さ */
  /* バレル厚 */
  rc = (rp + 0.001) + ac;

  /* 最下層コイル半径 */
  /*  Total_turn=0; */
  /* Num_lay=4;          %Number of layer */
  /* lpn=10;             %Number of lay per num */
  /* コイルに必要なデータ取り出し方 */
  /* コイル巻き数をNcとすると、iをNcまでふった時に、すべてのコイルにアクセスできる必要がある。 */
  /* つまり、コイルの指定は、単一のアドレスで引き出す必要がある */
  /* 必要なデータは、半径データ、ワイヤ径、z座標である。 */
  /* z座標は、zc_wからの相対位置が示されていればよい。 */
  /* 半径は、ワイヤ径を考慮して、矛盾なきよう。 */
  /* コイルデータをプロットして、矛盾がないか調べてみよう。 */
  i = obj_CoilMap->size[0] * obj_CoilMap->size[1];
  obj_CoilMap->size[0] = (int)(Num_lay * lpn);
  obj_CoilMap->size[1] = 3;
  emxEnsureCapacity_real_T(obj_CoilMap, i);
  loop_ub = (int)(Num_lay * lpn) * 3;
  for (i = 0; i < loop_ub; i++) {
    obj_CoilMap->data[i] = 0.0;
  }

  for (i = 0; i < (int)((Num_lay - 1.0) + 1.0); i++) {
    for (loop_ub = 0; loop_ub < (int)((lpn - 1.0) + 1.0); loop_ub++) {
      obj_CoilMap->data[(int)(((double)i * lpn + (double)loop_ub) + 1.0) - 1] =
        rc + ac * (double)i;
      obj_CoilMap->data[((int)(((double)i * lpn + (double)loop_ub) + 1.0) +
                         obj_CoilMap->size[0]) - 1] = ac * (double)loop_ub;

      /* この値にzc_wを加算すると、ワールド座標になる。 */
      obj_CoilMap->data[((int)(((double)i * lpn + (double)loop_ub) + 1.0) +
                         (obj_CoilMap->size[0] << 1)) - 1] = ac;
    }
  }

  Total = 0.0;
  for (i = 0; i < obj_CoilMap->size[0]; i++) {
    rc = obj_CoilMap->data[i + (obj_CoilMap->size[0] << 1)] * 0.5;

    /* Cupper wire resistance */
    /*      temp=20; %温度ケルビン */
    /*      length=5.9; %長さ　メートル */
    /*      dia=0.5*1e-3; %半径 メートル */
    /* Cu_roh=(1)*1.68e-8; */
    Total += 1.8097632000000003E-8 * (obj_CoilMap->data[i] * 2.0 *
      3.1415926535897931) / (rc * rc * 3.1415926535897931);
  }

  rc = 0.0;
  for (i = 0; i < obj_CoilMap->size[0]; i++) {
    for (loop_ub = 0; loop_ub < obj_CoilMap->size[0]; loop_ub++) {
      if (1 + i != 1 + loop_ub) {
        rc += MutualInd(obj_CoilMap->data[i], obj_CoilMap->data[loop_ub], fabs
                        (obj_CoilMap->data[i + obj_CoilMap->size[0]] -
                         obj_CoilMap->data[loop_ub + obj_CoilMap->size[0]]));
      } else {
        /* Self induntance */
        /* Umoto method */
        /* meter */
        /* 変数への入力はすべてメートル */
        rc += 1.2566370614359173E-6 * obj_CoilMap->data[i] / 4.0 +
          1.2566370614359173E-6 * obj_CoilMap->data[i] * (log(8.0 *
          obj_CoilMap->data[i] / (obj_CoilMap->data[i + (obj_CoilMap->size[0] <<
          1)] * 0.5)) - 2.0);
      }
    }
  }

  /* obj.Zc=zc; */
  *obj_NumOfTurn = lpn * Num_lay;
  *obj_SelfInd = rc;
  *obj_Resist = Total;
  *obj_Zc = 0.0;
}

static void __anon_fcn(double CP_Lsw, double CP_LD, double Coil_NumOfTurn, const
  emxArray_real_T *Coil_CoilMap, double Proj_NumOfElem, const emxArray_real_T
  *Proj_ElemMap, double Proj_Mass, const emxArray_real_T *AA, const
  emxArray_real_T *bb, const emxArray_real_T *y, emxArray_real_T *varargout_1)
{
  emxArray_real_T *b_Proj_ElemMap;
  int i;
  int loop_ub;
  emxArray_real_T *b_AA;
  emxArray_real_T *b_bb;
  emxArray_real_T *Mpc_vec;
  double K;
  int m;
  emxArray_int32_T *r2;
  emxArray_real_T *b_y;
  emxArray_real_T *A;
  int aoffset;
  emxInit_real_T(&b_Proj_ElemMap, 2);
  i = b_Proj_ElemMap->size[0] * b_Proj_ElemMap->size[1];
  b_Proj_ElemMap->size[0] = Proj_ElemMap->size[0];
  b_Proj_ElemMap->size[1] = 3;
  emxEnsureCapacity_real_T(b_Proj_ElemMap, i);
  loop_ub = Proj_ElemMap->size[0] * Proj_ElemMap->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_Proj_ElemMap->data[i] = Proj_ElemMap->data[i];
  }

  emxInit_real_T(&b_AA, 2);
  i = b_AA->size[0] * b_AA->size[1];
  b_AA->size[0] = AA->size[0];
  b_AA->size[1] = AA->size[1];
  emxEnsureCapacity_real_T(b_AA, i);
  loop_ub = AA->size[0] * AA->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_AA->data[i] = AA->data[i];
  }

  emxInit_real_T(&b_bb, 2);
  i = b_bb->size[0] * b_bb->size[1];
  b_bb->size[0] = bb->size[0];
  b_bb->size[1] = bb->size[1];
  emxEnsureCapacity_real_T(b_bb, i);
  loop_ub = bb->size[0] * bb->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_bb->data[i] = bb->data[i];
  }

  emxInit_real_T1(&Mpc_vec, 1);

  /* dydt = zeros(3,1); */
  K = y->data[(int)((2.0 + Proj_NumOfElem) + 3.0) - 1];

  /* v=0; */
  /* AA(2,2)=CP.Lsw+CP.LD; */
  /* OdeFuncに書く内容　試験的にこっちに書いてる　ここ、速度が要求されるので、高速な実行を得たい。 */
  /* Mutual_p2cをベクトルで返すやつ　逐一呼ぶより、効率がよかろうと思ってな。 */
  i = Mpc_vec->size[0];
  Mpc_vec->size[0] = (int)Proj_NumOfElem;
  emxEnsureCapacity_real_T1(Mpc_vec, i);
  loop_ub = (int)Proj_NumOfElem;
  for (i = 0; i < loop_ub; i++) {
    Mpc_vec->data[i] = 0.0;
  }

  for (loop_ub = 0; loop_ub < (int)Proj_NumOfElem; loop_ub++) {
    for (i = 0; i < (int)Coil_NumOfTurn; i++) {
      Mpc_vec->data[loop_ub] += MutualInd(Coil_CoilMap->data[i],
        b_Proj_ElemMap->data[loop_ub], (Coil_CoilMap->data[i +
        Coil_CoilMap->size[0]] - b_Proj_ElemMap->data[loop_ub +
        b_Proj_ElemMap->size[0]]) - K);
    }
  }

  if (3.0 > 2.0 + Proj_NumOfElem) {
    i = 0;
    m = 0;
  } else {
    i = 2;
    m = (int)(2.0 + Proj_NumOfElem);
  }

  emxInit_int32_T(&r2, 1);
  loop_ub = r2->size[0];
  r2->size[0] = m - i;
  emxEnsureCapacity_int32_T(r2, loop_ub);
  loop_ub = m - i;
  for (m = 0; m < loop_ub; m++) {
    r2->data[m] = i + m;
  }

  loop_ub = r2->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_AA->data[b_AA->size[0] * r2->data[i]] = Mpc_vec->data[i];
  }

  if (3.0 > 2.0 + Proj_NumOfElem) {
    i = 0;
  } else {
    i = 2;
  }

  loop_ub = Mpc_vec->size[0];
  for (m = 0; m < loop_ub; m++) {
    b_AA->data[i + m] = Mpc_vec->data[m];
  }

  /* FF=Mutual_p2c_vec(Coil,Proj); */
  /* 試験的に書いてる */
  /* Coil.Zc=6.9e-3; */
  /* OdeFuncに書く内容　試験的にこっちに書いてる　ここ、速度が要求されるので、高速な実行を得たい。 */
  /* これ、関数ハンドルうまく使って、値をうまく渡したい。一回だけ実行して変数に書き込む技を身に着けたい。もうグローバル変数でも良い気がするが。 */
  /* vを掛ける */
  Mutual_dx_p2c_vec(Coil_NumOfTurn, Coil_CoilMap, Proj_NumOfElem, b_Proj_ElemMap,
                    K, Mpc_vec);
  emxFree_real_T(&b_Proj_ElemMap);
  if (3.0 > 2.0 + Proj_NumOfElem) {
    i = 0;
  } else {
    i = 2;
  }

  K = -y->data[(int)((2.0 + Proj_NumOfElem) + 2.0) - 1];
  loop_ub = Mpc_vec->size[0];
  for (m = 0; m < loop_ub; m++) {
    b_bb->data[i + m] = K * Mpc_vec->data[m];
  }

  if (3.0 > 2.0 + Proj_NumOfElem) {
    i = 0;
    m = 0;
  } else {
    i = 2;
    m = (int)(2.0 + Proj_NumOfElem);
  }

  loop_ub = r2->size[0];
  r2->size[0] = m - i;
  emxEnsureCapacity_int32_T(r2, loop_ub);
  loop_ub = m - i;
  for (m = 0; m < loop_ub; m++) {
    r2->data[m] = i + m;
  }

  emxInit_real_T1(&b_y, 1);
  K = -y->data[(int)((2.0 + Proj_NumOfElem) + 2.0) - 1];
  i = b_y->size[0];
  b_y->size[0] = Mpc_vec->size[0];
  emxEnsureCapacity_real_T1(b_y, i);
  loop_ub = Mpc_vec->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_y->data[i] = K * Mpc_vec->data[i];
  }

  loop_ub = r2->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_bb->data[b_bb->size[0] * r2->data[i]] = b_y->data[i];
  }

  K = y->data[0] / Proj_Mass;
  if (3.0 > 2.0 + Proj_NumOfElem) {
    i = 0;
    m = 0;
  } else {
    i = 2;
    m = (int)(2.0 + Proj_NumOfElem);
  }

  loop_ub = r2->size[0];
  r2->size[0] = m - i;
  emxEnsureCapacity_int32_T(r2, loop_ub);
  loop_ub = m - i;
  for (m = 0; m < loop_ub; m++) {
    r2->data[m] = i + m;
  }

  i = b_y->size[0];
  b_y->size[0] = Mpc_vec->size[0];
  emxEnsureCapacity_real_T1(b_y, i);
  loop_ub = Mpc_vec->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_y->data[i] = K * Mpc_vec->data[i];
  }

  loop_ub = r2->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_bb->data[((int)(4.0 + Proj_NumOfElem) + b_bb->size[0] * r2->data[i]) - 1] =
      b_y->data[i];
  }

  emxFree_real_T(&b_y);
  emxFree_int32_T(&r2);
  if (y->data[(int)((2.0 + Proj_NumOfElem) + 1.0) - 1] > 0.0) {
    b_AA->data[1 + b_AA->size[0]] = CP_Lsw + 100.0;

    /* bb(2+1+Proj.NumOfElem,1)=-1/CP.C; */
    /* bb(2+1+Proj.NumOfElem,2)=-1/CP.C; */
  } else {
    b_AA->data[1 + b_AA->size[0]] = CP_Lsw + CP_LD;

    /* bb(2+1+Proj.NumOfElem,1)=0; */
    /* bb(2+1+Proj.NumOfElem,2)=0; */
  }

  emxInit_real_T(&A, 2);
  if ((b_AA->size[0] == 0) || (b_AA->size[1] == 0)) {
    i = A->size[0] * A->size[1];
    A->size[0] = b_AA->size[0];
    A->size[1] = b_AA->size[1];
    emxEnsureCapacity_real_T(A, i);
    loop_ub = b_AA->size[0] * b_AA->size[1];
    for (i = 0; i < loop_ub; i++) {
      A->data[i] = b_AA->data[i];
    }
  } else {
    invNxN(b_AA, A);
  }

  emxFree_real_T(&b_AA);
  if ((b_bb->size[1] == 1) || (y->size[0] == 1)) {
    i = Mpc_vec->size[0];
    Mpc_vec->size[0] = b_bb->size[0];
    emxEnsureCapacity_real_T1(Mpc_vec, i);
    loop_ub = b_bb->size[0];
    for (i = 0; i < loop_ub; i++) {
      Mpc_vec->data[i] = 0.0;
      aoffset = b_bb->size[1];
      for (m = 0; m < aoffset; m++) {
        Mpc_vec->data[i] += b_bb->data[i + b_bb->size[0] * m] * y->data[m];
      }
    }
  } else {
    m = b_bb->size[0];
    i = Mpc_vec->size[0];
    Mpc_vec->size[0] = b_bb->size[0];
    emxEnsureCapacity_real_T1(Mpc_vec, i);
    for (i = 1; i <= m; i++) {
      Mpc_vec->data[i - 1] = 0.0;
    }

    for (loop_ub = 0; loop_ub < b_bb->size[1]; loop_ub++) {
      if (y->data[loop_ub] != 0.0) {
        aoffset = loop_ub * m;
        for (i = 0; i < m; i++) {
          Mpc_vec->data[i] += y->data[loop_ub] * b_bb->data[aoffset + i];
        }
      }
    }
  }

  emxFree_real_T(&b_bb);
  if ((A->size[1] == 1) || (Mpc_vec->size[0] == 1)) {
    i = varargout_1->size[0];
    varargout_1->size[0] = A->size[0];
    emxEnsureCapacity_real_T1(varargout_1, i);
    loop_ub = A->size[0];
    for (i = 0; i < loop_ub; i++) {
      varargout_1->data[i] = 0.0;
      aoffset = A->size[1];
      for (m = 0; m < aoffset; m++) {
        varargout_1->data[i] += A->data[i + A->size[0] * m] * Mpc_vec->data[m];
      }
    }
  } else {
    m = A->size[0];
    i = varargout_1->size[0];
    varargout_1->size[0] = A->size[0];
    emxEnsureCapacity_real_T1(varargout_1, i);
    for (i = 1; i <= m; i++) {
      varargout_1->data[i - 1] = 0.0;
    }

    for (loop_ub = 0; loop_ub < A->size[1]; loop_ub++) {
      if (Mpc_vec->data[loop_ub] != 0.0) {
        aoffset = loop_ub * m;
        for (i = 0; i < m; i++) {
          varargout_1->data[i] += Mpc_vec->data[loop_ub] * A->data[aoffset + i];
        }
      }
    }
  }

  emxFree_real_T(&A);
  emxFree_real_T(&Mpc_vec);
}

static void callODEFunctionNSM(double c_odeFcn_tunableEnvironment_f1_, double
  odeFcn_tunableEnvironment_f1_LD, double c_odeFcn_tunableEnvironment_f2_, const
  emxArray_real_T *d_odeFcn_tunableEnvironment_f2_, double
  c_odeFcn_tunableEnvironment_f3_, const emxArray_real_T
  *d_odeFcn_tunableEnvironment_f3_, double e_odeFcn_tunableEnvironment_f3_,
  const emxArray_real_T *odeFcn_tunableEnvironment_f4, const emxArray_real_T
  *odeFcn_tunableEnvironment_f5, const emxArray_real_T *y, emxArray_real_T *yp)
{
  __anon_fcn(c_odeFcn_tunableEnvironment_f1_, odeFcn_tunableEnvironment_f1_LD,
             c_odeFcn_tunableEnvironment_f2_, d_odeFcn_tunableEnvironment_f2_,
             c_odeFcn_tunableEnvironment_f3_, d_odeFcn_tunableEnvironment_f3_,
             e_odeFcn_tunableEnvironment_f3_, odeFcn_tunableEnvironment_f4,
             odeFcn_tunableEnvironment_f5, y, yp);
}

static double eps(double x)
{
  double r;
  double absxk;
  int exponent;
  absxk = fabs(x);
  if ((!rtIsInf(absxk)) && (!rtIsNaN(absxk))) {
    if (absxk <= 2.2250738585072014E-308) {
      r = 4.94065645841247E-324;
    } else {
      frexp(absxk, &exponent);
      r = ldexp(1.0, exponent - 53);
    }
  } else {
    r = rtNaN;
  }

  return r;
}

static void invNxN(const emxArray_real_T *x, emxArray_real_T *y)
{
  int n;
  int i1;
  int yk;
  emxArray_real_T *b_x;
  int b_n;
  emxArray_int32_T *ipiv;
  int k;
  int u1;
  emxArray_int32_T *p;
  int j;
  int mmj;
  int c;
  int ix;
  double smax;
  int jy;
  double s;
  int ijA;
  n = x->size[0];
  i1 = y->size[0] * y->size[1];
  y->size[0] = x->size[0];
  y->size[1] = x->size[1];
  emxEnsureCapacity_real_T(y, i1);
  yk = x->size[0] * x->size[1];
  for (i1 = 0; i1 < yk; i1++) {
    y->data[i1] = 0.0;
  }

  emxInit_real_T(&b_x, 2);
  i1 = b_x->size[0] * b_x->size[1];
  b_x->size[0] = x->size[0];
  b_x->size[1] = x->size[1];
  emxEnsureCapacity_real_T(b_x, i1);
  yk = x->size[0] * x->size[1];
  for (i1 = 0; i1 < yk; i1++) {
    b_x->data[i1] = x->data[i1];
  }

  yk = x->size[0];
  if (yk < 1) {
    b_n = 0;
  } else {
    b_n = yk;
  }

  emxInit_int32_T1(&ipiv, 2);
  i1 = ipiv->size[0] * ipiv->size[1];
  ipiv->size[0] = 1;
  ipiv->size[1] = b_n;
  emxEnsureCapacity_int32_T1(ipiv, i1);
  if (b_n > 0) {
    ipiv->data[0] = 1;
    yk = 1;
    for (k = 2; k <= b_n; k++) {
      yk++;
      ipiv->data[k - 1] = yk;
    }
  }

  if (x->size[0] < 1) {
    b_n = 0;
  } else {
    yk = x->size[0] - 1;
    u1 = x->size[0];
    if (yk < u1) {
      u1 = yk;
    }

    for (j = 0; j < u1; j++) {
      mmj = n - j;
      c = j * (n + 1);
      if (mmj < 1) {
        yk = -1;
      } else {
        yk = 0;
        if (mmj > 1) {
          ix = c;
          smax = fabs(b_x->data[c]);
          for (k = 2; k <= mmj; k++) {
            ix++;
            s = fabs(b_x->data[ix]);
            if (s > smax) {
              yk = k - 1;
              smax = s;
            }
          }
        }
      }

      if (b_x->data[c + yk] != 0.0) {
        if (yk != 0) {
          ipiv->data[j] = (j + yk) + 1;
          ix = j;
          yk += j;
          for (k = 1; k <= n; k++) {
            smax = b_x->data[ix];
            b_x->data[ix] = b_x->data[yk];
            b_x->data[yk] = smax;
            ix += n;
            yk += n;
          }
        }

        i1 = c + mmj;
        for (jy = c + 1; jy < i1; jy++) {
          b_x->data[jy] /= b_x->data[c];
        }
      }

      yk = n - j;
      b_n = (c + n) + 1;
      jy = c + n;
      for (k = 1; k < yk; k++) {
        smax = b_x->data[jy];
        if (b_x->data[jy] != 0.0) {
          ix = c + 1;
          i1 = mmj + b_n;
          for (ijA = b_n; ijA < i1 - 1; ijA++) {
            b_x->data[ijA] += b_x->data[ix] * -smax;
            ix++;
          }
        }

        jy += n;
        b_n += n;
      }
    }

    b_n = x->size[0];
  }

  emxInit_int32_T1(&p, 2);
  i1 = p->size[0] * p->size[1];
  p->size[0] = 1;
  p->size[1] = b_n;
  emxEnsureCapacity_int32_T1(p, i1);
  if (b_n > 0) {
    p->data[0] = 1;
    yk = 1;
    for (k = 2; k <= b_n; k++) {
      yk++;
      p->data[k - 1] = yk;
    }
  }

  for (k = 0; k < ipiv->size[1]; k++) {
    if (ipiv->data[k] > 1 + k) {
      yk = p->data[ipiv->data[k] - 1];
      p->data[ipiv->data[k] - 1] = p->data[k];
      p->data[k] = yk;
    }
  }

  emxFree_int32_T(&ipiv);
  for (k = 0; k < n; k++) {
    c = p->data[k] - 1;
    y->data[k + y->size[0] * (p->data[k] - 1)] = 1.0;
    for (j = k; j < n; j++) {
      if (y->data[j + y->size[0] * c] != 0.0) {
        for (jy = j + 1; jy < n; jy++) {
          y->data[jy + y->size[0] * c] -= y->data[j + y->size[0] * c] *
            b_x->data[jy + b_x->size[0] * j];
        }
      }
    }
  }

  emxFree_int32_T(&p);
  if ((x->size[0] == 0) || ((y->size[0] == 0) || (y->size[1] == 0))) {
  } else {
    for (j = 1; j <= n; j++) {
      yk = n * (j - 1) - 1;
      for (k = n; k > 0; k--) {
        b_n = n * (k - 1) - 1;
        if (y->data[k + yk] != 0.0) {
          y->data[k + yk] /= b_x->data[k + b_n];
          for (jy = 1; jy < k; jy++) {
            y->data[jy + yk] -= y->data[k + yk] * b_x->data[jy + b_n];
          }
        }
      }
    }
  }

  emxFree_real_T(&b_x);
}

static void maxAbsThresh(const emxArray_real_T *a, const emxArray_real_T *b,
  emxArray_real_T *y)
{
  unsigned int unnamed_idx_0;
  int k;
  int loop_ub;
  double b_a;
  double b_b;
  unnamed_idx_0 = (unsigned int)a->size[0];
  k = y->size[0];
  y->size[0] = (int)unnamed_idx_0;
  emxEnsureCapacity_real_T1(y, k);
  loop_ub = (int)unnamed_idx_0;
  for (k = 0; k < loop_ub; k++) {
    y->data[k] = 0.0;
  }

  for (k = 0; k < a->size[0]; k++) {
    b_a = fabs(a->data[k]);
    b_b = fabs(b->data[k]);
    if ((b_a > b_b) || rtIsNaN(b_b)) {
      if (b_a > 0.1) {
        b_b = b_a;
      } else {
        b_b = 0.1;
      }
    } else {
      if (!(b_b > 0.1)) {
        b_b = 0.1;
      }
    }

    y->data[k] = b_b;
  }
}

static double norm(const emxArray_real_T *x)
{
  double y;
  int k;
  double absx;
  if (x->size[0] == 0) {
    y = 0.0;
  } else {
    y = 0.0;
    for (k = 0; k < x->size[0]; k++) {
      absx = fabs(x->data[k]);
      if (rtIsNaN(absx) || (absx > y)) {
        y = absx;
      }
    }
  }

  return y;
}

static void ntrp45(const double t[3], double t0, const emxArray_real_T *b_y0,
                   double h, const emxArray_real_T *f, emxArray_real_T *y)
{
  emxArray_real_T *fhBI1;
  int i;
  int aoffset;
  emxArray_real_T *fhBI2;
  double b_y[7];
  static const double b[7] = { -2.859375, 0.0, 4.0431266846361185, -3.90625,
    2.7939268867924527, -1.5714285714285714, 1.5 };

  int m;
  int k;
  emxArray_real_T *fhBI3;
  static const double b_b[7] = { 3.0833333333333335, 0.0, -6.2893081761006293,
    10.416666666666666, -6.8773584905660377, 3.6666666666666665, -4.0 };

  emxArray_real_T *fhBI4;
  static const double c_b[7] = { -1.1328125, 0.0, 2.6954177897574123, -5.859375,
    3.7610554245283021, -1.9642857142857142, 2.5 };

  double s;
  emxInit_real_T1(&fhBI1, 1);
  i = f->size[0];
  aoffset = fhBI1->size[0];
  fhBI1->size[0] = i;
  emxEnsureCapacity_real_T1(fhBI1, aoffset);
  for (aoffset = 0; aoffset < i; aoffset++) {
    fhBI1->data[aoffset] = f->data[aoffset] * h;
  }

  for (i = 0; i < 7; i++) {
    b_y[i] = h * b[i];
  }

  emxInit_real_T1(&fhBI2, 1);
  m = f->size[0];
  aoffset = fhBI2->size[0];
  fhBI2->size[0] = f->size[0];
  emxEnsureCapacity_real_T1(fhBI2, aoffset);
  for (i = 1; i <= m; i++) {
    fhBI2->data[i - 1] = 0.0;
  }

  for (k = 0; k < 7; k++) {
    if (b_y[k] != 0.0) {
      aoffset = k * m;
      for (i = 0; i < m; i++) {
        fhBI2->data[i] += b_y[k] * f->data[aoffset + i];
      }
    }
  }

  for (i = 0; i < 7; i++) {
    b_y[i] = h * b_b[i];
  }

  emxInit_real_T1(&fhBI3, 1);
  m = f->size[0];
  aoffset = fhBI3->size[0];
  fhBI3->size[0] = f->size[0];
  emxEnsureCapacity_real_T1(fhBI3, aoffset);
  for (i = 1; i <= m; i++) {
    fhBI3->data[i - 1] = 0.0;
  }

  for (k = 0; k < 7; k++) {
    if (b_y[k] != 0.0) {
      aoffset = k * m;
      for (i = 0; i < m; i++) {
        fhBI3->data[i] += b_y[k] * f->data[aoffset + i];
      }
    }
  }

  for (i = 0; i < 7; i++) {
    b_y[i] = h * c_b[i];
  }

  emxInit_real_T1(&fhBI4, 1);
  m = f->size[0];
  aoffset = fhBI4->size[0];
  fhBI4->size[0] = f->size[0];
  emxEnsureCapacity_real_T1(fhBI4, aoffset);
  for (i = 1; i <= m; i++) {
    fhBI4->data[i - 1] = 0.0;
  }

  for (k = 0; k < 7; k++) {
    if (b_y[k] != 0.0) {
      aoffset = k * m;
      for (i = 0; i < m; i++) {
        fhBI4->data[i] += b_y[k] * f->data[aoffset + i];
      }
    }
  }

  aoffset = y->size[0] * y->size[1];
  y->size[0] = b_y0->size[0];
  y->size[1] = 3;
  emxEnsureCapacity_real_T(y, aoffset);
  for (i = 0; i < 3; i++) {
    s = (t[i] - t0) / h;
    for (k = 0; k < b_y0->size[0]; k++) {
      y->data[k + y->size[0] * i] = (((fhBI4->data[k] * s + fhBI3->data[k]) * s
        + fhBI2->data[k]) * s + fhBI1->data[k]) * s + b_y0->data[k];
    }
  }

  emxFree_real_T(&fhBI4);
  emxFree_real_T(&fhBI3);
  emxFree_real_T(&fhBI2);
  emxFree_real_T(&fhBI1);
}

static void ode45(double ode_tunableEnvironment_f1_Lsw, double
                  ode_tunableEnvironment_f1_LD, double
                  c_ode_tunableEnvironment_f2_Num, const emxArray_real_T
                  *c_ode_tunableEnvironment_f2_Coi, double
                  c_ode_tunableEnvironment_f3_Num, const emxArray_real_T
                  *c_ode_tunableEnvironment_f3_Ele, double
                  ode_tunableEnvironment_f3_Mass, const emxArray_real_T
                  *ode_tunableEnvironment_f4, const emxArray_real_T
                  *ode_tunableEnvironment_f5, const double tspan[2], const
                  emxArray_real_T *b_y0, emxArray_real_T *varargout_1,
                  emxArray_real_T *varargout_2)
{
  emxArray_real_T *f0;
  double tfinal;
  int neq;
  int c_y0[1];
  emxArray_real_T d_y0;
  int y;
  emxArray_real_T *tout;
  int i0;
  emxArray_real_T *yout;
  int aoffset;
  int nout;
  int ix;
  double twidth;
  double safehmax;
  double d0;
  double err;
  double hmax;
  double hmin;
  double absh;
  emxArray_real_T *b_y;
  unsigned int y0_idx_0;
  int iy;
  emxArray_real_T *maxval;
  unsigned int y_idx_0;
  emxArray_real_T *r0;
  double t;
  emxArray_real_T *f;
  double tdir;
  boolean_T MinStepExit;
  boolean_T Done;
  emxArray_real_T *varargin_1;
  emxArray_real_T *result;
  cell_wrap_4 reshapes[2];
  emxArray_real_T *r1;
  emxArray_real_T *b_yout;
  int exitg1;
  boolean_T NoFailedAttempts;
  int exitg2;
  int Bcolidx;
  int j;
  double tnew;
  int m;
  static const double x[21] = { 0.2, 0.075, 0.225, 0.97777777777777775,
    -3.7333333333333334, 3.5555555555555554, 2.9525986892242035,
    -11.595793324188385, 9.8228928516994358, -0.29080932784636487,
    2.8462752525252526, -10.757575757575758, 8.9064227177434727,
    0.27840909090909088, -0.2735313036020583, 0.091145833333333329, 0.0,
    0.44923629829290207, 0.65104166666666663, -0.322376179245283,
    0.13095238095238096 };

  int outidx;
  static const double B[7] = { 0.0012326388888888888, 0.0,
    -0.0042527702905061394, 0.036979166666666667, -0.05086379716981132,
    0.0419047619047619, -0.025 };

  double toutnew[4];
  double tref[3];
  boolean_T empty_non_axis_sizes;
  double tmp_data[200];
  emxInit_real_T1(&f0, 1);
  tfinal = tspan[1];
  neq = b_y0->size[1];
  c_y0[0] = b_y0->size[1];
  d_y0 = *b_y0;
  d_y0.size = (int *)&c_y0;
  d_y0.numDimensions = 1;
  callODEFunctionNSM(ode_tunableEnvironment_f1_Lsw, ode_tunableEnvironment_f1_LD,
                     c_ode_tunableEnvironment_f2_Num,
                     c_ode_tunableEnvironment_f2_Coi,
                     c_ode_tunableEnvironment_f3_Num,
                     c_ode_tunableEnvironment_f3_Ele,
                     ode_tunableEnvironment_f3_Mass, ode_tunableEnvironment_f4,
                     ode_tunableEnvironment_f5, &d_y0, f0);
  y = (int)(8192.0 / (double)b_y0->size[1]) + 4;
  if (200 < y) {
    y = 200;
  }

  emxInit_real_T(&tout, 2);
  i0 = tout->size[0] * tout->size[1];
  tout->size[0] = 1;
  tout->size[1] = y;
  emxEnsureCapacity_real_T(tout, i0);
  for (i0 = 0; i0 < y; i0++) {
    tout->data[i0] = 0.0;
  }

  emxInit_real_T(&yout, 2);
  i0 = yout->size[0] * yout->size[1];
  yout->size[0] = b_y0->size[1];
  yout->size[1] = y;
  emxEnsureCapacity_real_T(yout, i0);
  aoffset = b_y0->size[1] * y;
  for (i0 = 0; i0 < aoffset; i0++) {
    yout->data[i0] = 0.0;
  }

  nout = 1;
  tout->data[0] = 0.0;
  ix = b_y0->size[1];
  for (i0 = 0; i0 < ix; i0++) {
    yout->data[i0] = b_y0->data[i0];
  }

  twidth = fabs(tspan[1]);
  safehmax = fabs(tspan[1]);
  if (rtIsNaN(safehmax)) {
    d0 = 0.0;
  } else {
    d0 = safehmax;
  }

  safehmax = 3.5527136788005009E-15 * d0;
  err = 0.1 * twidth;
  if ((err > safehmax) || rtIsNaN(safehmax)) {
    safehmax = err;
  }

  if ((twidth < safehmax) || rtIsNaN(safehmax)) {
    hmax = twidth;
  } else {
    hmax = safehmax;
  }

  hmin = 16.0 * eps(0.0);
  safehmax = fabs(tspan[1]);
  if ((hmax < safehmax) || rtIsNaN(safehmax)) {
    absh = hmax;
  } else {
    absh = safehmax;
  }

  emxInit_real_T1(&b_y, 1);
  y0_idx_0 = (unsigned int)b_y0->size[1];
  i0 = b_y->size[0];
  b_y->size[0] = (int)y0_idx_0;
  emxEnsureCapacity_real_T1(b_y, i0);
  for (iy = 0; iy < b_y0->size[1]; iy++) {
    b_y->data[iy] = fabs(b_y0->data[iy]);
  }

  emxInit_real_T1(&maxval, 1);
  y0_idx_0 = (unsigned int)b_y->size[0];
  y_idx_0 = (unsigned int)b_y->size[0];
  i0 = maxval->size[0];
  maxval->size[0] = (int)y_idx_0;
  emxEnsureCapacity_real_T1(maxval, i0);
  for (iy = 0; iy < (int)y0_idx_0; iy++) {
    err = b_y->data[iy];
    if (!(err > 0.1)) {
      err = 0.1;
    }

    maxval->data[iy] = err;
  }

  emxInit_real_T1(&r0, 1);
  rdivide(f0, maxval, r0);
  safehmax = norm(r0) / 0.12679145539688907;
  if (absh * safehmax > 1.0) {
    absh = 1.0 / safehmax;
  }

  if (!((absh > hmin) || rtIsNaN(hmin))) {
    absh = hmin;
  }

  t = 0.0;
  i0 = b_y->size[0];
  b_y->size[0] = b_y0->size[1];
  emxEnsureCapacity_real_T1(b_y, i0);
  aoffset = b_y0->size[1];
  for (i0 = 0; i0 < aoffset; i0++) {
    b_y->data[i0] = b_y0->data[i0];
  }

  emxInit_real_T(&f, 2);
  i0 = f->size[0] * f->size[1];
  f->size[0] = b_y0->size[1];
  f->size[1] = 7;
  emxEnsureCapacity_real_T(f, i0);
  aoffset = b_y0->size[1] * 7;
  for (i0 = 0; i0 < aoffset; i0++) {
    f->data[i0] = 0.0;
  }

  aoffset = f0->size[0];
  for (i0 = 0; i0 < aoffset; i0++) {
    f->data[i0] = f0->data[i0];
  }

  tdir = tspan[1];
  if (tspan[1] < 0.0) {
    tdir = -1.0;
  } else if (tspan[1] > 0.0) {
    tdir = 1.0;
  } else {
    if (tspan[1] == 0.0) {
      tdir = 0.0;
    }
  }

  MinStepExit = false;
  Done = false;
  emxInit_real_T(&varargin_1, 2);
  emxInit_real_T(&result, 2);
  emxInitMatrix_cell_wrap_4(reshapes);
  emxInit_real_T1(&r1, 1);
  emxInit_real_T(&b_yout, 2);
  do {
    exitg1 = 0;
    hmin = 16.0 * eps(t);
    if ((hmin > absh) || rtIsNaN(absh)) {
      safehmax = hmin;
    } else {
      safehmax = absh;
    }

    if ((hmax < safehmax) || rtIsNaN(safehmax)) {
      absh = hmax;
    } else {
      absh = safehmax;
    }

    safehmax = tdir * absh;
    if (1.1 * absh >= fabs(tfinal - t)) {
      safehmax = tfinal - t;
      absh = fabs(safehmax);
      Done = true;
    }

    NoFailedAttempts = true;
    do {
      exitg2 = 0;
      Bcolidx = 5;
      for (j = 0; j < 5; j++) {
        Bcolidx += j;
        i0 = f0->size[0];
        f0->size[0] = b_y->size[0];
        emxEnsureCapacity_real_T1(f0, i0);
        aoffset = b_y->size[0];
        for (i0 = 0; i0 < aoffset; i0++) {
          f0->data[i0] = b_y->data[i0];
        }

        if ((neq != 0) && (!(safehmax == 0.0))) {
          ix = Bcolidx + 1;
          aoffset = neq * j;
          for (m = 1; m <= aoffset + 1; m += neq) {
            twidth = safehmax * x[ix - 6];
            iy = 0;
            i0 = (m + neq) - 1;
            for (outidx = m; outidx <= i0; outidx++) {
              f0->data[iy] += f->data[outidx - 1] * twidth;
              iy++;
            }

            ix++;
          }
        }

        callODEFunctionNSM(ode_tunableEnvironment_f1_Lsw,
                           ode_tunableEnvironment_f1_LD,
                           c_ode_tunableEnvironment_f2_Num,
                           c_ode_tunableEnvironment_f2_Coi,
                           c_ode_tunableEnvironment_f3_Num,
                           c_ode_tunableEnvironment_f3_Ele,
                           ode_tunableEnvironment_f3_Mass,
                           ode_tunableEnvironment_f4, ode_tunableEnvironment_f5,
                           f0, r0);
        aoffset = r0->size[0];
        for (i0 = 0; i0 < aoffset; i0++) {
          f->data[i0 + f->size[0] * (j + 1)] = r0->data[i0];
        }
      }

      tnew = t + safehmax;
      if (Done) {
        tnew = tfinal;
      }

      i0 = f0->size[0];
      f0->size[0] = b_y->size[0];
      emxEnsureCapacity_real_T1(f0, i0);
      aoffset = b_y->size[0];
      for (i0 = 0; i0 < aoffset; i0++) {
        f0->data[i0] = b_y->data[i0];
      }

      if ((neq != 0) && (!(safehmax == 0.0))) {
        aoffset = neq * 5;
        for (m = 1; m <= aoffset + 1; m += neq) {
          twidth = safehmax * x[Bcolidx];
          iy = 0;
          i0 = (m + neq) - 1;
          for (outidx = m; outidx <= i0; outidx++) {
            f0->data[iy] += f->data[outidx - 1] * twidth;
            iy++;
          }

          Bcolidx++;
        }
      }

      callODEFunctionNSM(ode_tunableEnvironment_f1_Lsw,
                         ode_tunableEnvironment_f1_LD,
                         c_ode_tunableEnvironment_f2_Num,
                         c_ode_tunableEnvironment_f2_Coi,
                         c_ode_tunableEnvironment_f3_Num,
                         c_ode_tunableEnvironment_f3_Ele,
                         ode_tunableEnvironment_f3_Mass,
                         ode_tunableEnvironment_f4, ode_tunableEnvironment_f5,
                         f0, r0);
      aoffset = r0->size[0];
      for (i0 = 0; i0 < aoffset; i0++) {
        f->data[i0 + f->size[0] * 6] = r0->data[i0];
      }

      m = f->size[0];
      i0 = maxval->size[0];
      maxval->size[0] = f->size[0];
      emxEnsureCapacity_real_T1(maxval, i0);
      for (ix = 1; ix <= m; ix++) {
        maxval->data[ix - 1] = 0.0;
      }

      for (iy = 0; iy < 7; iy++) {
        if (B[iy] != 0.0) {
          aoffset = iy * m;
          for (ix = 0; ix < m; ix++) {
            maxval->data[ix] += B[iy] * f->data[aoffset + ix];
          }
        }
      }

      maxAbsThresh(b_y, f0, r0);
      rdivide(maxval, r0, r1);
      err = absh * norm(r1);
      if (err > 0.0001) {
        if (absh <= hmin) {
          MinStepExit = true;
          exitg2 = 1;
        } else {
          if (NoFailedAttempts) {
            NoFailedAttempts = false;
            safehmax = 0.8 * rt_powd_snf(0.0001 / err, 0.2);
            if (0.1 > safehmax) {
              safehmax = 0.1;
            }

            safehmax *= absh;
            if ((hmin > safehmax) || rtIsNaN(safehmax)) {
              absh = hmin;
            } else {
              absh = safehmax;
            }
          } else {
            safehmax = 0.5 * absh;
            if ((hmin > safehmax) || rtIsNaN(safehmax)) {
              absh = hmin;
            } else {
              absh = safehmax;
            }
          }

          safehmax = tdir * absh;
          Done = false;
        }
      } else {
        exitg2 = 1;
      }
    } while (exitg2 == 0);

    if (MinStepExit) {
      exitg1 = 1;
    } else {
      outidx = nout;
      safehmax = tnew - t;
      for (i0 = 0; i0 < 3; i0++) {
        twidth = t + safehmax * (0.25 + 0.25 * (double)i0);
        toutnew[i0] = twidth;
        tref[i0] = twidth;
      }

      toutnew[3] = tnew;
      ntrp45(tref, t, b_y, tnew - t, f, varargin_1);
      if (!(varargin_1->size[0] == 0)) {
        m = varargin_1->size[0];
      } else if (!(f0->size[0] == 0)) {
        m = f0->size[0];
      } else {
        m = 0;
      }

      empty_non_axis_sizes = (m == 0);
      if (empty_non_axis_sizes || (!(varargin_1->size[0] == 0))) {
        aoffset = 3;
      } else {
        aoffset = 0;
      }

      if (empty_non_axis_sizes || (!(f0->size[0] == 0))) {
        ix = 1;
      } else {
        ix = 0;
      }

      i0 = result->size[0] * result->size[1];
      result->size[0] = m;
      result->size[1] = aoffset + ix;
      emxEnsureCapacity_real_T(result, i0);
      for (i0 = 0; i0 < aoffset; i0++) {
        for (iy = 0; iy < m; iy++) {
          result->data[iy + result->size[0] * i0] = varargin_1->data[iy + m * i0];
        }
      }

      for (i0 = 0; i0 < ix; i0++) {
        for (iy = 0; iy < m; iy++) {
          result->data[iy + result->size[0] * (i0 + aoffset)] = f0->data[iy + m *
            i0];
        }
      }

      nout += 4;
      if (nout > tout->size[1]) {
        ix = tout->size[1];
        i0 = tout->size[0] * tout->size[1];
        tout->size[1] = ix + y;
        emxEnsureCapacity_real_T(tout, i0);
        if (0 <= y - 1) {
          memset(&tmp_data[0], 0, (unsigned int)(y * (int)sizeof(double)));
        }

        for (i0 = 0; i0 < y; i0++) {
          tout->data[ix + i0] = tmp_data[i0];
        }

        if (!((yout->size[0] == 0) || (yout->size[1] == 0))) {
          ix = yout->size[0];
        } else if (!((neq == 0) || (y == 0))) {
          ix = neq;
        } else {
          ix = yout->size[0];
          if (!(ix > 0)) {
            ix = 0;
          }

          if (neq > ix) {
            ix = neq;
          }
        }

        empty_non_axis_sizes = (ix == 0);
        if (empty_non_axis_sizes || (!((yout->size[0] == 0) || (yout->size[1] ==
               0)))) {
          m = yout->size[1];
        } else {
          m = 0;
        }

        if (empty_non_axis_sizes || (!((neq == 0) || (y == 0)))) {
          aoffset = y;
        } else {
          aoffset = 0;
        }

        i0 = reshapes[1].f1->size[0] * reshapes[1].f1->size[1];
        reshapes[1].f1->size[0] = ix;
        reshapes[1].f1->size[1] = aoffset;
        emxEnsureCapacity_real_T(reshapes[1].f1, i0);
        aoffset *= ix;
        for (i0 = 0; i0 < aoffset; i0++) {
          reshapes[1].f1->data[i0] = 0.0;
        }

        i0 = b_yout->size[0] * b_yout->size[1];
        b_yout->size[0] = ix;
        b_yout->size[1] = m + reshapes[1].f1->size[1];
        emxEnsureCapacity_real_T(b_yout, i0);
        for (i0 = 0; i0 < m; i0++) {
          for (iy = 0; iy < ix; iy++) {
            b_yout->data[iy + b_yout->size[0] * i0] = yout->data[iy + ix * i0];
          }
        }

        aoffset = reshapes[1].f1->size[1];
        for (i0 = 0; i0 < aoffset; i0++) {
          ix = reshapes[1].f1->size[0];
          for (iy = 0; iy < ix; iy++) {
            b_yout->data[iy + b_yout->size[0] * (i0 + m)] = reshapes[1].f1->
              data[iy + reshapes[1].f1->size[0] * i0];
          }
        }

        i0 = yout->size[0] * yout->size[1];
        yout->size[0] = b_yout->size[0];
        yout->size[1] = b_yout->size[1];
        emxEnsureCapacity_real_T(yout, i0);
        aoffset = b_yout->size[1];
        for (i0 = 0; i0 < aoffset; i0++) {
          ix = b_yout->size[0];
          for (iy = 0; iy < ix; iy++) {
            yout->data[iy + yout->size[0] * i0] = b_yout->data[iy + b_yout->
              size[0] * i0];
          }
        }
      }

      for (iy = 0; iy < 4; iy++) {
        tout->data[iy + outidx] = toutnew[iy];
        for (j = 0; j < neq; j++) {
          yout->data[j + yout->size[0] * (iy + outidx)] = result->data[j +
            result->size[0] * iy];
        }
      }

      if (Done) {
        exitg1 = 1;
      } else {
        if (NoFailedAttempts) {
          safehmax = 1.25 * rt_powd_snf(err / 0.0001, 0.2);
          if (safehmax > 0.2) {
            absh /= safehmax;
          } else {
            absh *= 5.0;
          }
        }

        t = tnew;
        i0 = b_y->size[0];
        b_y->size[0] = f0->size[0];
        emxEnsureCapacity_real_T1(b_y, i0);
        aoffset = f0->size[0];
        for (i0 = 0; i0 < aoffset; i0++) {
          b_y->data[i0] = f0->data[i0];
        }

        ix = f->size[0] - 1;
        for (i0 = 0; i0 <= ix; i0++) {
          f->data[i0] = f->data[i0 + f->size[0] * 6];
        }
      }
    }
  } while (exitg1 == 0);

  emxFree_real_T(&b_yout);
  emxFree_real_T(&r1);
  emxFree_real_T(&r0);
  emxFreeMatrix_cell_wrap_4(reshapes);
  emxFree_real_T(&result);
  emxFree_real_T(&varargin_1);
  emxFree_real_T(&maxval);
  emxFree_real_T(&b_y);
  emxFree_real_T(&f);
  emxFree_real_T(&f0);
  if (1 > nout) {
    aoffset = 0;
  } else {
    aoffset = nout;
  }

  i0 = varargout_1->size[0];
  varargout_1->size[0] = aoffset;
  emxEnsureCapacity_real_T1(varargout_1, i0);
  for (i0 = 0; i0 < aoffset; i0++) {
    varargout_1->data[i0] = tout->data[i0];
  }

  emxFree_real_T(&tout);
  if (1 > nout) {
    aoffset = 0;
  } else {
    aoffset = nout;
  }

  ix = yout->size[0];
  i0 = varargout_2->size[0] * varargout_2->size[1];
  varargout_2->size[0] = aoffset;
  varargout_2->size[1] = ix;
  emxEnsureCapacity_real_T(varargout_2, i0);
  for (i0 = 0; i0 < ix; i0++) {
    for (iy = 0; iy < aoffset; iy++) {
      varargout_2->data[iy + varargout_2->size[0] * i0] = yout->data[i0 +
        yout->size[0] * iy];
    }
  }

  emxFree_real_T(&yout);
}

static void rdivide(const emxArray_real_T *x, const emxArray_real_T *y,
                    emxArray_real_T *z)
{
  int i2;
  int loop_ub;
  i2 = z->size[0];
  z->size[0] = x->size[0];
  emxEnsureCapacity_real_T1(z, i2);
  loop_ub = x->size[0];
  for (i2 = 0; i2 < loop_ub; i2++) {
    z->data[i2] = x->data[i2] / y->data[i2];
  }
}

static double rt_powd_snf(double u0, double u1)
{
  double y;
  double d1;
  double d2;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d1 = fabs(u0);
    d2 = fabs(u1);
    if (rtIsInf(u1)) {
      if (d1 == 1.0) {
        y = 1.0;
      } else if (d1 > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d2 == 0.0) {
      y = 1.0;
    } else if (d2 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

double LIM_FreeE(double Rsw, double Lsw, double RD, double LD, double C, double
                 Rci, double Vc_init, double rp, double ac, double Num_layC,
                 double lpnC, double ap, double Num_layP, double lpnP, double zp,
                 double t_end)
{
  double Eff;
  emxArray_real_T *Coil_CoilMap;
  emxArray_real_T *Proj_ElemMap;
  emxArray_real_T *Proj_ResistMap;
  emxArray_real_T *Proj_SelfIndMap;
  double Coil_NumOfTurn;
  double Coil_SelfInd;
  double Coil_Resist;
  double Nc;
  int i;
  int loop_ub;
  double ProjMass;
  emxArray_real_T *AA;
  double rb;
  emxArray_real_T *bb;
  emxArray_real_T *b_y0;
  emxArray_real_T *tunableEnvironment_f2_CoilMap;
  emxArray_real_T *tunableEnvironment_f3_ElemMap;
  emxArray_real_T *tunableEnvironment_f4;
  emxArray_real_T *tunableEnvironment_f5;
  emxArray_real_T *c_this_tunableEnvironment_f2_Co;
  emxArray_real_T *c_this_tunableEnvironment_f3_El;
  emxArray_real_T *this_tunableEnvironment_f4;
  emxArray_real_T *this_tunableEnvironment_f5;
  double dv0[2];
  emxInit_real_T(&Coil_CoilMap, 2);
  emxInit_real_T(&Proj_ElemMap, 2);
  emxInit_real_T1(&Proj_ResistMap, 1);
  emxInit_real_T1(&Proj_SelfIndMap, 1);

  /* CoilGene機能に、コイルパラメータ配列を出力する機能を付加しましょう。 */
  RectCoil_RectCoil(rp, ac, Num_layC, lpnC, &Coil_NumOfTurn, &Coil_SelfInd,
                    &Coil_Resist, Coil_CoilMap, &Nc);

  /* Projectle data gene */
  /* Proj(最外半径、分解能、層数、層あたり巻き数) */
  /* 矩形断面プロジェクタイルデータ作成 */
  /* rp=20*1e-3;         %スリーブ最外径 */
  /* ap=1*1e-3;          %スリーブ分解能 */
  /* Num_lay=3;          %Number of layer */
  /* lpn=10;             %Number of lay per num */
  /* コイルに必要なデータ取り出し方と同じ */
  /* コイル巻き数をNcとすると、iをNcまでふった時に、すべてのコイルにアクセスできる必要がある。 */
  /* つまり、コイルの指定は、単一のアドレスで引き出す必要がある */
  /* 必要なデータは、半径データ、ワイヤ径、z座標である。 */
  /* z座標は、zc_wからの相対位置が示されていればよい。 */
  /* 半径は、ワイヤ径を考慮して、矛盾なきよう。 */
  /* コイルデータをプロットして、矛盾がないか調べてみよう。 */
  Nc = Num_layP * lpnP;
  i = Proj_ElemMap->size[0] * Proj_ElemMap->size[1];
  Proj_ElemMap->size[0] = (int)Nc;
  Proj_ElemMap->size[1] = 3;
  emxEnsureCapacity_real_T(Proj_ElemMap, i);
  loop_ub = (int)Nc * 3;
  for (i = 0; i < loop_ub; i++) {
    Proj_ElemMap->data[i] = 0.0;
  }

  for (i = 0; i < (int)((Num_layP - 1.0) + 1.0); i++) {
    for (loop_ub = 0; loop_ub < (int)((lpnP - 1.0) + 1.0); loop_ub++) {
      Proj_ElemMap->data[(int)(((double)i * lpnP + (double)loop_ub) + 1.0) - 1] =
        rp - ap * (double)i;

      /* 最外径を基準に、内側に積層する感じ */
      Proj_ElemMap->data[((int)(((double)i * lpnP + (double)loop_ub) + 1.0) +
                          Proj_ElemMap->size[0]) - 1] = ap * (double)loop_ub;

      /* この値にzp_wを加算すると、ワールド座標になる。 */
      Proj_ElemMap->data[((int)(((double)i * lpnP + (double)loop_ub) + 1.0) +
                          (Proj_ElemMap->size[0] << 1)) - 1] = ap;

      /* ワイヤ太さ（直径 */
      /*         Proj(i*lpn+j+1,4)=R_Al(Proj(i*lpn+j+1,1)^2*pi,ap,ap,20);         %R_Al(length,H,W,temp) */
      /*         Proj(i*lpn+j+1,5)=SelfInd(Proj(i*lpn+j+1,1),Proj(i*lpn+j+1,3)*0.5);  %SelfInd(環半径,ワイヤ半径); */
    }
  }

  /* Projectle data gene */
  i = Proj_SelfIndMap->size[0];
  Proj_SelfIndMap->size[0] = Proj_ElemMap->size[0];
  emxEnsureCapacity_real_T1(Proj_SelfIndMap, i);
  loop_ub = Proj_ElemMap->size[0];
  for (i = 0; i < loop_ub; i++) {
    Proj_SelfIndMap->data[i] = 0.0;
  }

  for (i = 0; i < Proj_ElemMap->size[0]; i++) {
    /* Self induntance */
    /* Umoto method */
    /* meter */
    /* 変数への入力はすべてメートル */
    Proj_SelfIndMap->data[i] += 1.2566370614359173E-6 * Proj_ElemMap->data[i] /
      4.0 + 1.2566370614359173E-6 * Proj_ElemMap->data[i] * (log(8.0 *
      Proj_ElemMap->data[i] / (Proj_ElemMap->data[i + (Proj_ElemMap->size[0] <<
      1)] * 0.5)) - 2.0);

    /* SelfInd(R,a) */
  }

  /* Projectle data gene */
  i = Proj_ResistMap->size[0];
  Proj_ResistMap->size[0] = Proj_ElemMap->size[0];
  emxEnsureCapacity_real_T1(Proj_ResistMap, i);
  loop_ub = Proj_ElemMap->size[0];
  for (i = 0; i < loop_ub; i++) {
    Proj_ResistMap->data[i] = 0.0;
  }

  for (i = 0; i < Proj_ElemMap->size[0]; i++) {
    /* R_Al */
    /* temp=10+273;　%温度ケルビン */
    /* length=1000;　%長さ　メートル */
    /* H  高さ */
    /* W　幅 */
    Proj_ResistMap->data[i] = 3.0399600000000006E-8 * (Proj_ElemMap->data[i] *
      2.0 * 3.1415926535897931) / (Proj_ElemMap->data[i + (Proj_ElemMap->size[0]
      << 1)] * Proj_ElemMap->data[i + (Proj_ElemMap->size[0] << 1)]);

    /* SelfInd(環半径,ワイヤ半径); */
  }

  /* R_Al(Proj(i*lpn+j+1,1)^2*pi,ap,ap,20); */
  /* R_Al(length,H,W,temp) */
  /* Projectle data gene */
  ProjMass = 0.0;
  for (i = 0; i < Proj_ElemMap->size[0]; i++) {
    Nc = Proj_ElemMap->data[i] + Proj_ElemMap->data[i + (Proj_ElemMap->size[0] <<
      1)] * 0.5;
    rb = Proj_ElemMap->data[i] - Proj_ElemMap->data[i + (Proj_ElemMap->size[0] <<
      1)] * 0.5;

    /* CalcMass */
    /* 外形、内径から円筒の断面積を計算、→掃引して体積算出→密度から質量変換 */
    /*  ra=50*1e-3; */
    /*  rb=0*1e-3; */
    /*  h=10e-3; */
    /* m^3→cm^3 */
    /* キログラム */
    ProjMass += 2.7 * ((Nc * Nc * 3.1415926535897931 - rb * rb *
                        3.1415926535897931) * Proj_ElemMap->data[i +
                       (Proj_ElemMap->size[0] << 1)] * 1.0E+6) * 0.001;

    /* SelfInd(環半径,ワイヤ半径); */
  }

  emxInit_real_T(&AA, 2);
  Nc = lpnP * Num_layP;

  /* obj.Zc=zc; */
  /* 位置マップと抵抗マップ、自己インダクタンスマップは別々に配置しよう。 */
  i = AA->size[0] * AA->size[1];
  AA->size[0] = (int)((Nc + 2.0) + 3.0);
  AA->size[1] = (int)((Nc + 2.0) + 3.0);
  emxEnsureCapacity_real_T(AA, i);
  loop_ub = (int)((Nc + 2.0) + 3.0) * (int)((Nc + 2.0) + 3.0);
  for (i = 0; i < loop_ub; i++) {
    AA->data[i] = 0.0;
  }

  emxInit_real_T(&bb, 2);
  i = bb->size[0] * bb->size[1];
  bb->size[0] = (int)((Nc + 2.0) + 3.0);
  bb->size[1] = (int)((Nc + 2.0) + 3.0);
  emxEnsureCapacity_real_T(bb, i);
  loop_ub = (int)((Nc + 2.0) + 3.0) * (int)((Nc + 2.0) + 3.0);
  for (i = 0; i < loop_ub; i++) {
    bb->data[i] = 0.0;
  }

  /* A行列の定数の入力 */
  for (i = 2; i - 2 < (int)((Nc + 2.0) + -2.0); i++) {
    AA->data[i + AA->size[0] * i] = Proj_SelfIndMap->data[(int)((3.0 + (double)
      (i - 2)) - 2.0) - 1];
    bb->data[i + bb->size[0] * i] = -Proj_ResistMap->data[(int)((3.0 + (double)
      (i - 2)) - 2.0) - 1];
  }

  emxFree_real_T(&Proj_SelfIndMap);
  AA->data[0] = Lsw + Coil_SelfInd;
  AA->data[AA->size[0]] = Lsw;
  AA->data[1] = Lsw;
  AA->data[1 + AA->size[0]] = Lsw + LD;
  AA->data[((int)(3.0 + Nc) + AA->size[0] * ((int)(3.0 + Nc) - 1)) - 1] = 1.0;
  AA->data[((int)(4.0 + Nc) + AA->size[0] * ((int)(4.0 + Nc) - 1)) - 1] = 1.0;
  AA->data[((int)(5.0 + Nc) + AA->size[0] * ((int)(5.0 + Nc) - 1)) - 1] = 1.0;

  /* これ、定数じゃね？→定数ですわこれ */
  for (i = 0; i < (int)Nc; i++) {
    for (loop_ub = 0; loop_ub < (int)Nc; loop_ub++) {
      if (1.0 + (double)i != 1.0 + (double)loop_ub) {
        AA->data[((int)(2.0 + (1.0 + (double)i)) + AA->size[0] * ((int)(2.0 +
                    (1.0 + (double)loop_ub)) - 1)) - 1] = MutualInd
          (Proj_ElemMap->data[i], Proj_ElemMap->data[loop_ub],
           Proj_ElemMap->data[i + Proj_ElemMap->size[0]] - Proj_ElemMap->
           data[loop_ub + Proj_ElemMap->size[0]]);
      }
    }
  }

  emxInit_real_T(&b_y0, 2);
  bb->data[0] = -((Rsw + Rci) + Coil_Resist);
  bb->data[bb->size[0]] = -Rsw;
  bb->data[1] = -Rsw;
  bb->data[1 + bb->size[0]] = -Rsw - RD;
  bb->data[(int)(3.0 + Nc) - 1] = -1.0 / C;
  bb->data[((int)(3.0 + Nc) + bb->size[0]) - 1] = -1.0 / C;
  bb->data[bb->size[0] * ((int)(3.0 + Nc) - 1)] = 1.0;
  bb->data[1 + bb->size[0] * ((int)(3.0 + Nc) - 1)] = 1.0;
  bb->data[((int)(5.0 + Nc) + bb->size[0] * ((int)(4.0 + Nc) - 1)) - 1] = 1.0;

  /* 試験的に書いてる */
  i = b_y0->size[0] * b_y0->size[1];
  b_y0->size[0] = 1;
  b_y0->size[1] = (int)((Nc + 2.0) + 3.0);
  emxEnsureCapacity_real_T(b_y0, i);
  loop_ub = (int)((Nc + 2.0) + 3.0);
  for (i = 0; i < loop_ub; i++) {
    b_y0->data[i] = 0.0;
  }

  emxInit_real_T(&tunableEnvironment_f2_CoilMap, 2);
  b_y0->data[0] = 1.0;
  b_y0->data[(int)((Nc + 2.0) + 3.0) - 1] = zp;
  b_y0->data[(int)((Nc + 2.0) + 1.0) - 1] = Vc_init;

  /* PCplot */
  i = tunableEnvironment_f2_CoilMap->size[0] *
    tunableEnvironment_f2_CoilMap->size[1];
  tunableEnvironment_f2_CoilMap->size[0] = Coil_CoilMap->size[0];
  tunableEnvironment_f2_CoilMap->size[1] = 3;
  emxEnsureCapacity_real_T(tunableEnvironment_f2_CoilMap, i);
  loop_ub = Coil_CoilMap->size[0] * Coil_CoilMap->size[1];
  for (i = 0; i < loop_ub; i++) {
    tunableEnvironment_f2_CoilMap->data[i] = Coil_CoilMap->data[i];
  }

  emxFree_real_T(&Coil_CoilMap);
  emxInit_real_T(&tunableEnvironment_f3_ElemMap, 2);
  i = tunableEnvironment_f3_ElemMap->size[0] *
    tunableEnvironment_f3_ElemMap->size[1];
  tunableEnvironment_f3_ElemMap->size[0] = Proj_ElemMap->size[0];
  tunableEnvironment_f3_ElemMap->size[1] = 3;
  emxEnsureCapacity_real_T(tunableEnvironment_f3_ElemMap, i);
  loop_ub = Proj_ElemMap->size[0] * Proj_ElemMap->size[1];
  for (i = 0; i < loop_ub; i++) {
    tunableEnvironment_f3_ElemMap->data[i] = Proj_ElemMap->data[i];
  }

  emxFree_real_T(&Proj_ElemMap);
  emxInit_real_T(&tunableEnvironment_f4, 2);
  i = tunableEnvironment_f4->size[0] * tunableEnvironment_f4->size[1];
  tunableEnvironment_f4->size[0] = AA->size[0];
  tunableEnvironment_f4->size[1] = AA->size[1];
  emxEnsureCapacity_real_T(tunableEnvironment_f4, i);
  loop_ub = AA->size[0] * AA->size[1];
  for (i = 0; i < loop_ub; i++) {
    tunableEnvironment_f4->data[i] = AA->data[i];
  }

  emxInit_real_T(&tunableEnvironment_f5, 2);
  i = tunableEnvironment_f5->size[0] * tunableEnvironment_f5->size[1];
  tunableEnvironment_f5->size[0] = bb->size[0];
  tunableEnvironment_f5->size[1] = bb->size[1];
  emxEnsureCapacity_real_T(tunableEnvironment_f5, i);
  loop_ub = bb->size[0] * bb->size[1];
  for (i = 0; i < loop_ub; i++) {
    tunableEnvironment_f5->data[i] = bb->data[i];
  }

  emxFree_real_T(&bb);
  emxInit_real_T(&c_this_tunableEnvironment_f2_Co, 2);
  i = c_this_tunableEnvironment_f2_Co->size[0] *
    c_this_tunableEnvironment_f2_Co->size[1];
  c_this_tunableEnvironment_f2_Co->size[0] = tunableEnvironment_f2_CoilMap->
    size[0];
  c_this_tunableEnvironment_f2_Co->size[1] = 3;
  emxEnsureCapacity_real_T(c_this_tunableEnvironment_f2_Co, i);
  loop_ub = tunableEnvironment_f2_CoilMap->size[0] *
    tunableEnvironment_f2_CoilMap->size[1];
  for (i = 0; i < loop_ub; i++) {
    c_this_tunableEnvironment_f2_Co->data[i] =
      tunableEnvironment_f2_CoilMap->data[i];
  }

  emxFree_real_T(&tunableEnvironment_f2_CoilMap);
  emxInit_real_T(&c_this_tunableEnvironment_f3_El, 2);
  i = c_this_tunableEnvironment_f3_El->size[0] *
    c_this_tunableEnvironment_f3_El->size[1];
  c_this_tunableEnvironment_f3_El->size[0] = tunableEnvironment_f3_ElemMap->
    size[0];
  c_this_tunableEnvironment_f3_El->size[1] = 3;
  emxEnsureCapacity_real_T(c_this_tunableEnvironment_f3_El, i);
  loop_ub = tunableEnvironment_f3_ElemMap->size[0] *
    tunableEnvironment_f3_ElemMap->size[1];
  for (i = 0; i < loop_ub; i++) {
    c_this_tunableEnvironment_f3_El->data[i] =
      tunableEnvironment_f3_ElemMap->data[i];
  }

  emxFree_real_T(&tunableEnvironment_f3_ElemMap);
  emxInit_real_T(&this_tunableEnvironment_f4, 2);
  i = this_tunableEnvironment_f4->size[0] * this_tunableEnvironment_f4->size[1];
  this_tunableEnvironment_f4->size[0] = tunableEnvironment_f4->size[0];
  this_tunableEnvironment_f4->size[1] = tunableEnvironment_f4->size[1];
  emxEnsureCapacity_real_T(this_tunableEnvironment_f4, i);
  loop_ub = tunableEnvironment_f4->size[0] * tunableEnvironment_f4->size[1];
  for (i = 0; i < loop_ub; i++) {
    this_tunableEnvironment_f4->data[i] = tunableEnvironment_f4->data[i];
  }

  emxFree_real_T(&tunableEnvironment_f4);
  emxInit_real_T(&this_tunableEnvironment_f5, 2);
  i = this_tunableEnvironment_f5->size[0] * this_tunableEnvironment_f5->size[1];
  this_tunableEnvironment_f5->size[0] = tunableEnvironment_f5->size[0];
  this_tunableEnvironment_f5->size[1] = tunableEnvironment_f5->size[1];
  emxEnsureCapacity_real_T(this_tunableEnvironment_f5, i);
  loop_ub = tunableEnvironment_f5->size[0] * tunableEnvironment_f5->size[1];
  for (i = 0; i < loop_ub; i++) {
    this_tunableEnvironment_f5->data[i] = tunableEnvironment_f5->data[i];
  }

  emxFree_real_T(&tunableEnvironment_f5);
  dv0[0] = 0.0;
  dv0[1] = t_end;
  ode45(Lsw, LD, Coil_NumOfTurn, c_this_tunableEnvironment_f2_Co, Nc,
        c_this_tunableEnvironment_f3_El, ProjMass, this_tunableEnvironment_f4,
        this_tunableEnvironment_f5, dv0, b_y0, Proj_ResistMap, AA);

  /* plot(t,y(:,Proj.NumOfElem+2+1)) */
  /* ploter; */
  Eff = 100.0 * (0.5 * ProjMass * (AA->data[(AA->size[0] + AA->size[0] * ((int)
    ((2.0 + Nc) + 2.0) - 1)) - 1] * AA->data[(AA->size[0] + AA->size[0] * ((int)
    ((2.0 + Nc) + 2.0) - 1)) - 1])) / (0.5 * C * (Vc_init * Vc_init));

  /*   */
  emxFree_real_T(&this_tunableEnvironment_f5);
  emxFree_real_T(&this_tunableEnvironment_f4);
  emxFree_real_T(&c_this_tunableEnvironment_f3_El);
  emxFree_real_T(&c_this_tunableEnvironment_f2_Co);
  emxFree_real_T(&b_y0);
  emxFree_real_T(&AA);
  emxFree_real_T(&Proj_ResistMap);
  return Eff;
}

void LIM_FreeE_initialize(void)
{
  rt_InitInfAndNaN(8U);
}

void LIM_FreeE_terminate(void)
{
  /* (no terminate code required) */
}

/* End of code generation (LIM_FreeE.c) */
