#include <Novice.h>
#define _USE_MATH_DEFINES
#include "../../../KamataEngine/DirectXGame/math/Vector4.h"
#include <imgui.h>
#include <stdint.h>

#include "AffineMatrix.h"

const char kWindowTitle[] = "LE2B_27_ヤマダ_ナオ_3_2_確認課題";

const float kWindowWidth = 1280.0f;
const float kWindowHeight = 720.0f;

/// <summary>
/// 行列の和
/// </summary>
/// <param name="m1"></param>
/// <param name="m2"></param>
/// <returns></returns>
Matrix4x4 Add(const Matrix4x4 &m1, const Matrix4x4 &m2) {
  Matrix4x4 result;

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      result.m[i][j] = m1.m[i][j] + m2.m[i][j];
    }
  }

  return result;
}

/// <summary>
/// 行列の差
/// </summary>
/// <param name="m1"></param>
/// <param name="m2"></param>
/// <returns></returns>
Matrix4x4 Subtract(const Matrix4x4 &m1, const Matrix4x4 &m2) {
  Matrix4x4 result;

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      result.m[i][j] = m1.m[i][j] - m2.m[i][j];
    }
  }

  return result;
}

/// <summary>
/// 3次元ベクトルの和
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
Vector3 Add(const Vector3 &v1, const Vector3 v2) {
  Vector3 result;
  result.x = v1.x + v2.x;
  result.y = v1.y + v2.y;
  result.z = v1.z + v2.z;

  return result;
}

/// <summary>
/// 3次元ベクトルの差
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
Vector3 Subtract(const Vector3 v1, const Vector3 v2) {
  Vector3 result;
  result.x = v1.x - v2.x;
  result.y = v1.y - v2.y;
  result.z = v1.z - v2.z;

  return result;
}

/// <summary>
/// 3次元ベクトルの積
/// </summary>
/// <param name="s"></param>
/// <param name="v"></param>
/// <returns></returns>
Vector3 Multiply(float s, const Vector3 &v) {
  Vector3 result;
  result.x = v.x * s;
  result.y = v.y * s;
  result.z = v.z * s;

  return result;
}

/// <summary>
/// X軸周りの回転行列
/// </summary>
/// <param name="radian"></param>
/// <returns></returns>
Matrix4x4 MakeRotateXMatrix(float radian) {
  Matrix4x4 result;
  result.m[0][0] = 1.0f;
  result.m[0][1] = 0.0f;
  result.m[0][2] = 0.0f;
  result.m[0][3] = 0.0f;
  result.m[1][0] = 0.0f;
  result.m[1][1] = std::cos(radian);
  result.m[1][2] = std::sin(radian);
  result.m[1][3] = 0.0f;
  result.m[2][0] = 0.0f;
  result.m[2][1] = -std::sin(radian);
  result.m[2][2] = std::cos(radian);
  result.m[2][3] = 0.0f;
  result.m[3][0] = 0.0f;
  result.m[3][1] = 0.0f;
  result.m[3][2] = 0.0f;
  result.m[3][3] = 1.0f;

  return result;
}

/// <summary>
/// Y軸周りの回転行列
/// </summary>
/// <param name="radian"></param>
/// <returns></returns>
Matrix4x4 MakeRotateYMatrix(float radian) {
  Matrix4x4 result;
  result.m[0][0] = std::cos(radian);
  result.m[0][1] = 0.0f;
  result.m[0][2] = -std::sin(radian);
  result.m[0][3] = 0.0f;
  result.m[1][0] = 0.0f;
  result.m[1][1] = 1.0f;
  result.m[1][2] = 0.0f;
  result.m[1][3] = 0.0f;
  result.m[2][0] = std::sin(radian);
  result.m[2][1] = 0.0f;
  result.m[2][2] = std::cos(radian);
  result.m[2][3] = 0.0f;
  result.m[3][0] = 0.0f;
  result.m[3][1] = 0.0f;
  result.m[3][2] = 0.0f;
  result.m[3][3] = 1.0f;

  return result;
}

/// <summary>
/// Z軸周りの回転行列
/// </summary>
/// <param name="radian"></param>
/// <returns></returns>
Matrix4x4 MakeRotateZMatrix(float radian) {
  Matrix4x4 result;
  result.m[0][0] = std::cos(radian);
  result.m[0][1] = std::sin(radian);
  result.m[0][2] = 0.0f;
  result.m[0][3] = 0.0f;
  result.m[1][0] = -std::sin(radian);
  result.m[1][1] = std::cos(radian);
  result.m[1][2] = 0.0f;
  result.m[1][3] = 0.0f;
  result.m[2][0] = 0.0f;
  result.m[2][1] = 0.0f;
  result.m[2][2] = 1.0f;
  result.m[2][3] = 0.0f;
  result.m[3][0] = 0.0f;
  result.m[3][1] = 0.0f;
  result.m[3][2] = 0.0f;
  result.m[3][3] = 1.0f;

  return result;
}

/// ========================================
/// 演算子のオーバーロード
/// ======================================== 

/// <summary>
/// 3次元ベクトルの和
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
Vector3 operator+(const Vector3 &v1, const Vector3 &v2) { return Add(v1, v2); }

/// <summary>
/// 3次元ベクトルの差
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
Vector3 operator-(const Vector3 &v1, const Vector3 &v2) {
  return Subtract(v1, v2);
}

/// <summary>
/// 3次元ベクトルの積
/// </summary>
/// <param name="v"></param>
/// <param name="s"></param>
/// <returns></returns>
Vector3 operator*(const Vector3 &v, float s) { return Multiply(s, v); }

/// <summary>
/// 3次元ベクトルの商
/// </summary>
/// <param name="v"></param>
/// <param name="s"></param>
/// <returns></returns>
Vector3 operator/(const Vector3 &v, float s) { return Multiply(1.0f / s, v); }

/// <summary>
/// 行列の和
/// </summary>
/// <param name="m1"></param>
/// <param name="m2"></param>
/// <returns></returns>
Matrix4x4 operator+(const Matrix4x4 &m1, const Matrix4x4 &m2) {
  return Add(m1, m2);
}

/// <summary>
/// 行列の差
/// </summary>
/// <param name="m1"></param>
/// <param name="m2"></param>
/// <returns></returns>
Matrix4x4 operator-(const Matrix4x4 &m1, const Matrix4x4 &m2) {
  return Subtract(m1, m2);
}

/// <summary>
/// 行列の積
/// </summary>
/// <param name="m1"></param>
/// <param name="m2"></param>
/// <returns></returns>
Matrix4x4 operator*(const Matrix4x4 m1, const Matrix4x4 &m2) {
  return Multiply(m1, m2);
}

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

  // ライブラリの初期化
  Novice::Initialize(kWindowTitle, static_cast<int>(kWindowWidth),
                     static_cast<int>(kWindowHeight));

  // キー入力結果を受け取る箱
  char keys[256] = {0};
  char preKeys[256] = {0};

  Vector3 a{0.2f, 1.0f, 0.0f};
  Vector3 b{2.4f, 3.1f, 1.2f};
  Vector3 c = a + b;
  Vector3 d = a - b;
  Vector3 e = a * 2.4f;

  Vector3 rotate{0.4f, 1.43f, -0.8f};
  Matrix4x4 rotateXMatrix = MakeRotateXMatrix(rotate.x);
  Matrix4x4 rotateYMatrix = MakeRotateYMatrix(rotate.y);
  Matrix4x4 rotateZMatrix = MakeRotateZMatrix(rotate.z);
  Matrix4x4 rotateMatrix = rotateXMatrix * rotateYMatrix * rotateZMatrix;

  // ウィンドウの×ボタンが押されるまでループ
  while (Novice::ProcessMessage() == 0) {
    // フレームの開始
    Novice::BeginFrame();

    // キー入力を受け取る
    memcpy(preKeys, keys, 256);
    Novice::GetHitKeyStateAll(keys);

    ///
    /// ↓更新処理ここから
    ///

    ImGui::Begin("Window");
    ImGui::Text("c: %f, %f %f", c.x, c.y, c.z);
    ImGui::Text("d: %f, %f %f", d.x, d.y, d.z);
    ImGui::Text("e: %f, %f %f", e.x, e.y, e.z);
    ImGui::Text(
        "matrix:\n%f, %f, %f, %f\n%f, %f, %f, %f\n%f, %f, %f, %f\n%f, %f, %f, %f",
        rotateMatrix.m[0][0], rotateMatrix.m[0][1], rotateMatrix.m[0][2],
        rotateMatrix.m[0][3], rotateMatrix.m[1][0], rotateMatrix.m[1][1],
        rotateMatrix.m[1][2], rotateMatrix.m[1][3], rotateMatrix.m[2][0],
        rotateMatrix.m[2][1], rotateMatrix.m[2][2], rotateMatrix.m[2][3],
        rotateMatrix.m[3][0], rotateMatrix.m[3][1], rotateMatrix.m[3][2],
        rotateMatrix.m[3][3]);

    ImGui::End();

    ///
    /// ↑更新処理ここまで
    ///

    ///
    /// ↓描画処理ここから
    ///

    ///
    /// ↑描画処理ここまで
    ///

    // フレームの終了
    Novice::EndFrame();

    // ESCキーが押されたらループを抜ける
    if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
      break;
    }
  }

  // ライブラリの終了
  Novice::Finalize();
  return 0;
}