#include <Novice.h>
#define _USE_MATH_DEFINES
#include "../../../KamataEngine/DirectXGame/math/Vector4.h"
#include <imgui.h>
#include <stdint.h>

#include "AffineMatrix.h"

const char kWindowTitle[] = "LE2B_27_ヤマダ_ナオ_4_0_確認課題";

struct Spring {
  Vector3 anchor;           // 固定された端
  float naturalLength;      // 自然長
  float stiffness;          // ばね定数
  float dampingCoefficient; // 減衰係数
};

struct Ball {
  Vector3 position;     // 位置
  Vector3 velocity;     // 速度
  Vector3 acceleration; // 加速度
  float mass;           // 質量
  float radius;         // 半径
  unsigned int color;   // 色
};

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
/// 3次元ベクトルの積
/// </summary>
/// <param name="v"></param>
/// <param name="s"></param>
/// <returns></returns>
Vector3 operator*(float s, const Vector3 &v) { return v * s; }

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

/// <summary>
/// 長さの計算
/// </summary>
/// <param name="v"></param>
/// <returns></returns>
float Length(const Vector3 &v) {
  return std::sqrtf(std::powf(v.x, 2) + std::powf(v.y, 2) + std::powf(v.z, 2));
}

struct Transformation {
  Vector3 scale;
  Vector3 rotate;
  Vector3 translate;
};

/// <summary>
/// 逆行列
/// </summary>
/// <param name="m"></param>
/// <returns></returns>
Matrix4x4 Inverse(const Matrix4x4 &m) {
  Matrix4x4 result;
  float det;

  // 行列式を求めるための補助変数
  float a0 = m.m[0][0] * m.m[1][1] - m.m[0][1] * m.m[1][0];
  float a1 = m.m[0][0] * m.m[1][2] - m.m[0][2] * m.m[1][0];
  float a2 = m.m[0][0] * m.m[1][3] - m.m[0][3] * m.m[1][0];
  float a3 = m.m[0][1] * m.m[1][2] - m.m[0][2] * m.m[1][1];
  float a4 = m.m[0][1] * m.m[1][3] - m.m[0][3] * m.m[1][1];
  float a5 = m.m[0][2] * m.m[1][3] - m.m[0][3] * m.m[1][2];

  float b0 = m.m[2][0] * m.m[3][1] - m.m[2][1] * m.m[3][0];
  float b1 = m.m[2][0] * m.m[3][2] - m.m[2][2] * m.m[3][0];
  float b2 = m.m[2][0] * m.m[3][3] - m.m[2][3] * m.m[3][0];
  float b3 = m.m[2][1] * m.m[3][2] - m.m[2][2] * m.m[3][1];
  float b4 = m.m[2][1] * m.m[3][3] - m.m[2][3] * m.m[3][1];
  float b5 = m.m[2][2] * m.m[3][3] - m.m[2][3] * m.m[3][2];

  // 行列式
  det = a0 * b5 - a1 * b4 + a2 * b3 + a3 * b2 - a4 * b1 + a5 * b0;

  if (det == 0.0f) {
    // ゼロ行列を返す
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        result.m[i][j] = 0.0f;
    return result;
  }

  float invDet = 1.0f / det;

  // 逆行列を計算
  result.m[0][0] = (m.m[1][1] * b5 - m.m[1][2] * b4 + m.m[1][3] * b3) * invDet;
  result.m[0][1] = (-m.m[0][1] * b5 + m.m[0][2] * b4 - m.m[0][3] * b3) * invDet;
  result.m[0][2] = (m.m[3][1] * a5 - m.m[3][2] * a4 + m.m[3][3] * a3) * invDet;
  result.m[0][3] = (-m.m[2][1] * a5 + m.m[2][2] * a4 - m.m[2][3] * a3) * invDet;

  result.m[1][0] = (-m.m[1][0] * b5 + m.m[1][2] * b2 - m.m[1][3] * b1) * invDet;
  result.m[1][1] = (m.m[0][0] * b5 - m.m[0][2] * b2 + m.m[0][3] * b1) * invDet;
  result.m[1][2] = (-m.m[3][0] * a5 + m.m[3][2] * a2 - m.m[3][3] * a1) * invDet;
  result.m[1][3] = (m.m[2][0] * a5 - m.m[2][2] * a2 + m.m[2][3] * a1) * invDet;

  result.m[2][0] = (m.m[1][0] * b4 - m.m[1][1] * b2 + m.m[1][3] * b0) * invDet;
  result.m[2][1] = (-m.m[0][0] * b4 + m.m[0][1] * b2 - m.m[0][3] * b0) * invDet;
  result.m[2][2] = (m.m[3][0] * a4 - m.m[3][1] * a2 + m.m[3][3] * a0) * invDet;
  result.m[2][3] = (-m.m[2][0] * a4 + m.m[2][1] * a2 - m.m[2][3] * a0) * invDet;

  result.m[3][0] = (-m.m[1][0] * b3 + m.m[1][1] * b1 - m.m[1][2] * b0) * invDet;
  result.m[3][1] = (m.m[0][0] * b3 - m.m[0][1] * b1 + m.m[0][2] * b0) * invDet;
  result.m[3][2] = (-m.m[3][0] * a3 + m.m[3][1] * a1 - m.m[3][2] * a0) * invDet;
  result.m[3][3] = (m.m[2][0] * a3 - m.m[2][1] * a1 + m.m[2][2] * a0) * invDet;

  return result;
}

/// <summary>
/// 正射影行列(3次元版)
/// </summary>
/// <param name="left"></param>
/// <param name="top"></param>
/// <param name="right"></param>
/// <param name="bottom"></param>
/// <param name="nearClip"></param>
/// <param name="farClip"></param>
/// <returns></returns>
Matrix4x4 MakeOrthographicMatrix(float left, float top, float right,
                                 float bottom, float nearClip, float farClip) {
  Matrix4x4 result;
  result.m[0][0] = 2.0f / (right - left);
  result.m[0][1] = 0.0f;
  result.m[0][2] = 0.0f;
  result.m[0][3] = 0.0f;

  result.m[1][0] = 0.0f;
  result.m[1][1] = 2.0f / (top - bottom);
  result.m[1][2] = 0.0f;
  result.m[1][3] = 0.0f;

  result.m[2][0] = 0.0f;
  result.m[2][1] = 0.0f;
  result.m[2][2] = 1.0f / (farClip - nearClip);
  result.m[2][3] = 0.0f;

  result.m[3][0] = (left + right) / (left - right);
  result.m[3][1] = (top + bottom) / (bottom - top);
  result.m[3][2] = nearClip / (nearClip - farClip);
  result.m[3][3] = 1.0f;

  return result;
}

/// <summary>
/// 透視投影行列
/// </summary>
/// <param name="fovY"></param>
/// <param name="aspectRatio"></param>
/// <param name="nearClip"></param>
/// <param name="farClip"></param>
/// <returns></returns>
Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio,
                                   float nearClip, float farClip) {
  Matrix4x4 result;
  result.m[0][0] = 1.0f / aspectRatio * 1.0f / std::tan(fovY / 2.0f);
  result.m[0][1] = 0.0f;
  result.m[0][2] = 0.0f;
  result.m[0][3] = 0.0f;

  result.m[1][0] = 0.0f;
  result.m[1][1] = 1.0f / std::tan(fovY / 2.0f);
  result.m[1][2] = 0.0f;
  result.m[1][3] = 0.0f;

  result.m[2][0] = 0.0f;
  result.m[2][1] = 0.0f;
  result.m[2][2] = farClip / (farClip - nearClip);
  result.m[2][3] = 1.0f;

  result.m[3][0] = 0.0f;
  result.m[3][1] = 0.0f;
  result.m[3][2] = -nearClip * farClip / (farClip - nearClip);
  result.m[3][3] = 0.0f;

  return result;
}

/// <summary>
/// ビューポート変換行列
/// </summary>
/// <param name="left"></param>
/// <param name="top"></param>
/// <param name="width"></param>
/// <param name="height"></param>
/// <param name="minDepth"></param>
/// <param name="maxDepth"></param>
/// <returns></returns>
Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height,
                             float minDepth, float maxDepth) {
  Matrix4x4 result;
  result.m[0][0] = width / 2.0f;
  result.m[0][1] = 0.0f;
  result.m[0][2] = 0.0f;
  result.m[0][3] = 0.0f;

  result.m[1][0] = 0.0f;
  result.m[1][1] = -height / 2.0f;
  result.m[1][2] = 0.0f;
  result.m[1][3] = 0.0f;

  result.m[2][0] = 0.0f;
  result.m[2][1] = 0.0f;
  result.m[2][2] = maxDepth - minDepth;
  result.m[2][3] = 0.0f;

  result.m[3][0] = left + width / 2.0f;
  result.m[3][1] = top + height / 2.0f;
  result.m[3][2] = minDepth;
  result.m[3][3] = 1.0f;

  return result;
}

// 4成分ベクトル変換（W込み）
Vector4 Transform4(const Vector3 &v, const Matrix4x4 &m) {
  return {
      v.x * m.m[0][0] + v.y * m.m[1][0] + v.z * m.m[2][0] + 1.0f * m.m[3][0],
      v.x * m.m[0][1] + v.y * m.m[1][1] + v.z * m.m[2][1] + 1.0f * m.m[3][1],
      v.x * m.m[0][2] + v.y * m.m[1][2] + v.z * m.m[2][2] + 1.0f * m.m[3][2],
      v.x * m.m[0][3] + v.y * m.m[1][3] + v.z * m.m[2][3] + 1.0f * m.m[3][3]};
}

// NDC 空間への正規化 (Wで割る)
Vector3 ToNDC(const Vector4 &c) { return {c.x / c.w, c.y / c.w, c.z / c.w}; }

// 3成分ベクトル変換（ビューポート用）
Vector3 Transform(const Vector3 &v, const Matrix4x4 &m) {
  return {
      v.x * m.m[0][0] + v.y * m.m[1][0] + v.z * m.m[2][0] + 1.0f * m.m[3][0],
      v.x * m.m[0][1] + v.y * m.m[1][1] + v.z * m.m[2][1] + 1.0f * m.m[3][1],
      v.x * m.m[0][2] + v.y * m.m[1][2] + v.z * m.m[2][2] + 1.0f * m.m[3][2]};
}

/// <summary>
/// 内積
/// </summary>
/// <param name="a"></param>
/// <param name="b"></param>
/// <returns></returns>
float Dot(const Vector3 &a, const Vector3 &b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

/// <summary>
/// クロス積
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
Vector3 Cross(const Vector3 &v1, const Vector3 &v2) {
  return {v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
          v1.x * v2.y - v1.y * v2.x};
}

/// <summary>
/// カメラの操作
/// </summary>
/// <param name="camera"></param>
/// <param name="keys"></param>
/// <param name="wheel"></param>
/// <param name="mouseX"></param>
/// <param name="mouseY"></param>
/// <param name="prevMouseX"></param>
/// <param name="prevMouseY"></param>
void UpdateCameraControl(Transformation &camera, const char *keys, int wheel,
                         int mouseX, int mouseY, int &prevMouseX,
                         int &prevMouseY) {
  if (keys[DIK_D])
    camera.translate.x += 0.03f;
  if (keys[DIK_A])
    camera.translate.x -= 0.03f;
  if (keys[DIK_W])
    camera.translate.z += 0.03f;
  if (keys[DIK_S])
    camera.translate.z -= 0.03f;

  float deltaX = static_cast<float>(mouseX - prevMouseX);
  float deltaY = static_cast<float>(mouseY - prevMouseY);
  prevMouseX = mouseX;
  prevMouseY = mouseY;

  if (Novice::IsPressMouse(1)) {
    camera.rotate.y += deltaX * 0.005f;
    camera.rotate.x += deltaY * 0.005f;
    if (camera.rotate.x > 1.5f)
      camera.rotate.x = 1.5f;
    if (camera.rotate.x < -1.5f)
      camera.rotate.x = -1.5f;
  }

  camera.scale.z += static_cast<float>(wheel) * 0.001f;
  if (camera.scale.z < 0.120f)
    camera.scale.z = 0.120f;
}

Vector3 Normalize(const Vector3 &v) {
  Vector3 result;
  // ベクトルの長さを求める
  float length = Length(v);

  // ゼロ除算を避ける
  if (length != 0.0f) {
    result.x = v.x / length;
    result.y = v.y / length;
    result.z = v.z / length;
  } else {
    // 長さゼロならそのまま返す or 0ベクトル
    result = {0.0f, 0.0f, 0.0f};
  }
  return result;
}

/// <summary>
/// グリッド線を描画
/// </summary>
/// <param name="viewProjectionMatrix"></param>
/// <param name="viewportMatrix"></param>
void DrawGrid(const Matrix4x4 &viewProjectionMatrix,
              const Matrix4x4 &viewportMatrix) {
  const float kGridHalfWidth = 2.0f;
  const uint32_t kSubdivision = 10;
  const float kGridEvery =
      (kGridHalfWidth * 2.0f) / static_cast<float>(kSubdivision);

  for (uint32_t xIndex = 0; xIndex <= kSubdivision; ++xIndex) {
    float x = -kGridHalfWidth + xIndex * kGridEvery;

    Vector3 start = {x, 0.0f, -kGridHalfWidth};
    Vector3 end = {x, 0.0f, kGridHalfWidth};

    Vector4 ndcStart4 = Transform4(start, viewProjectionMatrix);
    Vector4 ndcEnd4 = Transform4(end, viewProjectionMatrix);

    if (ndcStart4.w != 0.0f && ndcEnd4.w != 0.0f) {
      Vector3 ndcStart = {ndcStart4.x / ndcStart4.w, ndcStart4.y / ndcStart4.w,
                          ndcStart4.z / ndcStart4.w};
      Vector3 ndcEnd = {ndcEnd4.x / ndcEnd4.w, ndcEnd4.y / ndcEnd4.w,
                        ndcEnd4.z / ndcEnd4.w};

      Vector3 screenStart = Transform(ndcStart, viewportMatrix);
      Vector3 screenEnd = Transform(ndcEnd, viewportMatrix);

      Novice::DrawLine(static_cast<int>(screenStart.x),
                       static_cast<int>(screenStart.y),
                       static_cast<int>(screenEnd.x),
                       static_cast<int>(screenEnd.y), 0xAAAAAAFF);
    }
  }

  for (uint32_t zIndex = 0; zIndex <= kSubdivision; ++zIndex) {
    float z = -kGridHalfWidth + zIndex * kGridEvery;

    Vector3 start = {-kGridHalfWidth, 0.0f, z};
    Vector3 end = {kGridHalfWidth, 0.0f, z};

    Vector4 ndcStart4 = Transform4(start, viewProjectionMatrix);
    Vector4 ndcEnd4 = Transform4(end, viewProjectionMatrix);

    if (ndcStart4.w != 0.0f && ndcEnd4.w != 0.0f) {
      Vector3 ndcStart = {ndcStart4.x / ndcStart4.w, ndcStart4.y / ndcStart4.w,
                          ndcStart4.z / ndcStart4.w};
      Vector3 ndcEnd = {ndcEnd4.x / ndcEnd4.w, ndcEnd4.y / ndcEnd4.w,
                        ndcEnd4.z / ndcEnd4.w};

      Vector3 screenStart = Transform(ndcStart, viewportMatrix);
      Vector3 screenEnd = Transform(ndcEnd, viewportMatrix);

      Novice::DrawLine(static_cast<int>(screenStart.x),
                       static_cast<int>(screenStart.y),
                       static_cast<int>(screenEnd.x),
                       static_cast<int>(screenEnd.y), 0xAAAAAAFF);
    }
  }
}

void DrawBall(const Ball &ball, const Matrix4x4 &viewProjectionMatrix,
              const Matrix4x4 &viewportMatrix) {
  // まずボールの3D位置を4Dに変換
  Vector4 ndcPos4 = Transform4(ball.position, viewProjectionMatrix);

  // W除算（NDC変換）
  if (ndcPos4.w == 0.0f) {
    return; // 非表示
  }
  Vector3 ndcPos = ToNDC(ndcPos4);

  // ビューポート変換
  Vector3 screenPos = Transform(ndcPos, viewportMatrix);

  // Novice::DrawEllipse は整数座標が必要なのでキャスト
  Novice::DrawEllipse(
      static_cast<int>(screenPos.x), static_cast<int>(screenPos.y),
      static_cast<int>(ball.radius), static_cast<int>(ball.radius), 0.0f,
      ball.color, kFillModeSolid);
}

void DrawSpringLine(const Spring &spring, const Ball &ball,
                    const Matrix4x4 &viewProjectionMatrix,
                    const Matrix4x4 &viewportMatrix) {
  // spring.anchor → 4D
  Vector4 ndcAnchor4 = Transform4(spring.anchor, viewProjectionMatrix);
  // W除算
  if (ndcAnchor4.w == 0.0f)
    return;
  Vector3 ndcAnchor = ToNDC(ndcAnchor4);
  // ビューポート
  Vector3 screenAnchor = Transform(ndcAnchor, viewportMatrix);

  // ball.position → 4D
  Vector4 ndcBall4 = Transform4(ball.position, viewProjectionMatrix);
  if (ndcBall4.w == 0.0f)
    return;
  Vector3 ndcBall = ToNDC(ndcBall4);
  Vector3 screenBall = Transform(ndcBall, viewportMatrix);

  // 描画
  Novice::DrawLine(
      static_cast<int>(screenAnchor.x), static_cast<int>(screenAnchor.y),
      static_cast<int>(screenBall.x), static_cast<int>(screenBall.y),
      0xFF0000FF // 赤
  );
}

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

  // ライブラリの初期化
  Novice::Initialize(kWindowTitle, static_cast<int>(kWindowWidth),
                     static_cast<int>(kWindowHeight));

  // キー入力結果を受け取る箱
  char keys[256] = {0};
  char preKeys[256] = {0};

  Spring spring = {
      {0.0f, 0.0f, 0.0f}, // 固定された端の位置
      1.0f,               // 自然長
      100.0f,             // ばね定数
      2.0f                // 減衰係数
  };

  Ball ball = {
      {1.2f, 0.0f, 0.0f}, // 位置
      {0.0f, 0.0f, 0.0f}, // 速度
      {0.0f, 0.0f, 0.0f}, // 加速度
      2.0f,               // 質量
      8.0f,               // 半径
      BLUE                // 色
  };

  Transformation camera{
      {1.0f, 1.0f, 1.0f},  // scale
      {0.26f, 0.0f, 0.0f}, // rotate
      {0.0f, 1.9f, -6.49f} // translate
  };

  int wheel = 0;
  int mouseX = 0;
  int mouseY = 0;
  int prevMouseX = 0;
  int prevMouseY = 0;

  Transformation transform{
      {1.0f, 1.0f, 1.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}};

  float deltaTime = 1.0f / 60.0f;

  int isMoving = false;

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

    if (ImGui::Button("start")) {
      if (!isMoving) {
        isMoving = true;
      }
    }

    ImGui::End();

    if (isMoving) {
      Vector3 diff = ball.position - spring.anchor;
      float length = Length(diff);
      if (length != 0.0f) {
        Vector3 direction = Normalize(diff);
        Vector3 restPosition = spring.anchor + direction * spring.naturalLength;
        Vector3 displacement = length * (ball.position - restPosition);
        Vector3 restoringForce = -spring.stiffness * displacement;

        // 減衰力
        Vector3 dampingForce = -spring.dampingCoefficient * ball.velocity;
        // 力が減衰する
        Vector3 force = restoringForce + dampingForce;

        ball.acceleration = force / ball.mass;
      }

      ball.velocity += ball.acceleration * deltaTime;
      ball.position += ball.velocity * deltaTime;
    }

    Novice::GetMousePosition(&mouseX, &mouseY);
    wheel = Novice::GetWheel();

    UpdateCameraControl(camera, keys, wheel, mouseX, mouseY, prevMouseX,
                        prevMouseY);

    Matrix4x4 worldMatrix = MakeAffineMatrix(transform.scale, transform.rotate,
                                             transform.translate);
    Matrix4x4 cameraMatrix =
        MakeAffineMatrix(camera.scale, camera.rotate, camera.translate);
    Matrix4x4 viewMatrix = Inverse(cameraMatrix);
    Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(
        0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
    Matrix4x4 viewProjectionMatrix = Multiply(viewMatrix, projectionMatrix);
    Matrix4x4 worldViewProjectionMatrix =
        Multiply(worldMatrix, viewProjectionMatrix);
    Matrix4x4 viewportMatrix = MakeViewportMatrix(
        0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);

    ///
    /// ↑更新処理ここまで
    ///

    ///
    /// ↓描画処理ここから
    ///

    DrawGrid(worldViewProjectionMatrix, viewportMatrix);

    // ばね
    DrawSpringLine(spring, ball, worldViewProjectionMatrix, viewportMatrix);

    // ばね先のボール
    DrawBall(ball, worldViewProjectionMatrix, viewportMatrix);

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