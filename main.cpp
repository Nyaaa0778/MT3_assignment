#include <Novice.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <imgui.h>
#include <stdint.h>
#include<algorithm>

const char kWindowTitle[] = "LE2B_27_ヤマダ_ナオ_2_4_確認課題";

const float kWindowWidth = 1280.0f;
const float kWindowHeight = 720.0f;

struct Vector3 {
  float x;
  float y;
  float z;
};

struct Vector4 {
  float x;
  float y;
  float z;
  float w;
};

struct Matrix4x4 {
  float m[4][4];
};

struct Camera {
  Vector3 scale;
  Vector3 rotate;
  Vector3 translate;
};

struct Sphere {
  Vector3 center;
  float radius;
};

struct AABB {
  Vector3 min;
  Vector3 max;
};

/// <summary>
/// ベクトルの加算
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
Vector3 Add(const Vector3 &v1, const Vector3 &v2) {
  Vector3 result;
  result.x = v1.x + v2.x;
  result.y = v1.y + v2.y;
  result.z = v1.z + v2.z;

  return result;
}

/// <summary>
/// 行列の積
/// </summary>
/// <param name="m1"></param>
/// <param name="m2"></param>
/// <returns></returns>
Matrix4x4 Multiply(const Matrix4x4 &m1, const Matrix4x4 &m2) {
  Matrix4x4 result;

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      result.m[i][j] = 0;
      for (int k = 0; k < 4; ++k) {
        result.m[i][j] += m1.m[i][k] * m2.m[k][j];
      }
    }
  }

  return result;
}

/// <summary>
/// 拡縮行列
/// </summary>
/// <param name="scale"></param>
/// <returns></returns>
Matrix4x4 MakeScaleMatrix(const Vector3 &scale) {
  Matrix4x4 result;

  result.m[0][0] = scale.x;
  result.m[0][1] = 0.0f;
  result.m[0][2] = 0.0f;
  result.m[0][3] = 0.0f;
  result.m[1][0] = 0.0f;
  result.m[1][1] = scale.y;
  result.m[1][2] = 0.0f;
  result.m[1][3] = 0.0f;
  result.m[2][0] = 0.0f;
  result.m[2][1] = 0.0f;
  result.m[2][2] = scale.z;
  result.m[2][3] = 0.0f;
  result.m[3][0] = 0.0f;
  result.m[3][1] = 0.0f;
  result.m[3][2] = 0.0f;
  result.m[3][3] = 1.0f;

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

/// <summary>
/// 平行移動行列
/// </summary>
/// <param name="translate"></param>
/// <returns></returns>
Matrix4x4 MakeTranslateMatrix(const Vector3 &translate) {
  Matrix4x4 result;
  result.m[0][0] = 1.0f;
  result.m[0][1] = 0.0f;
  result.m[0][2] = 0.0f;
  result.m[0][3] = 0.0f;
  result.m[1][0] = 0.0f;
  result.m[1][1] = 1.0f;
  result.m[1][2] = 0.0f;
  result.m[1][3] = 0.0f;
  result.m[2][0] = 0.0f;
  result.m[2][1] = 0.0f;
  result.m[2][2] = 1.0f;
  result.m[2][3] = 0.0f;
  result.m[3][0] = translate.x;
  result.m[3][1] = translate.y;
  result.m[3][2] = translate.z;
  result.m[3][3] = 1.0f;

  return result;
}

/// <summary>
/// アフィン変換行列
/// </summary>
/// <param name="scale">拡縮</param>
/// <param name="rotate">回転</param>
/// <param name="translate">移動</param>
/// <returns></returns>
Matrix4x4 MakeAffineMatrix(const Vector3 &scale, const Vector3 &rotate,
                           const Vector3 &translate) {

  Matrix4x4 scaleMatrix = MakeScaleMatrix(scale);
  Matrix4x4 rotateXMatrix = MakeRotateXMatrix(rotate.x);
  Matrix4x4 rotateYMatrix = MakeRotateYMatrix(rotate.y);
  Matrix4x4 rotateZMatrix = MakeRotateZMatrix(rotate.z);
  Matrix4x4 translateMatrix = MakeTranslateMatrix(translate);

  Matrix4x4 worldMatrix = Multiply(
      Multiply(Multiply(Multiply(scaleMatrix, rotateXMatrix), rotateYMatrix),
               rotateZMatrix),
      translateMatrix);

  return worldMatrix;
}

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
/// グリッド線の描画
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

void UpdateCameraControl(Camera &camera, const char *keys, int wheel,
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

float Length(const Vector3 &v1, const Vector3 &v2) {
  return sqrtf(powf(v1.x - v2.x, 2) + powf(v1.y - v2.y, 2) +
               powf(v1.z - v2.z, 2));
}

bool IsCollision(const AABB& aabb, const Sphere& sphere) {
  Vector3 closestPoint{std::clamp(sphere.center.x, aabb.min.x, aabb.max.x),
                       std::clamp(sphere.center.y, aabb.min.y, aabb.max.y),
                       std::clamp(sphere.center.z, aabb.min.z, aabb.max.z)};

  float distance = Length(closestPoint, sphere.center);

  return distance <= sphere.radius;
}

/// <summary>
/// 球体の描画
/// </summary>
/// <param name="sphere"></param>
/// <param name="viewProjectionMatrix"></param>
/// <param name="viewportMatrix"></param>
/// <param name="color"></param>
void DrawSphere(const Sphere &sphere, const Matrix4x4 &viewProjectionMatrix,
                const Matrix4x4 &viewportMatrix, uint32_t color) {
  const uint32_t kSubDivision = 16;
  const float kLonEvery = 2 * float(M_PI) / kSubDivision;
  const float kLatEvery = float(M_PI) / kSubDivision;

  for (uint32_t latIndex = 0; latIndex < kSubDivision; ++latIndex) {
    float lat = -float(M_PI) / 2.0f + kLatEvery * latIndex;
    float latNext = lat + kLatEvery;

    for (uint32_t lonIndex = 0; lonIndex < kSubDivision; ++lonIndex) {
      float lon = lonIndex * kLonEvery;
      float lonNext = lon + kLonEvery;

      Vector3 a = {
          sphere.center.x + sphere.radius * cosf(lat) * cosf(lon),
          sphere.center.y + sphere.radius * sinf(lat),
          sphere.center.z + sphere.radius * cosf(lat) * sinf(lon),
      };

      Vector3 b = {
          sphere.center.x + sphere.radius * cosf(latNext) * cosf(lon),
          sphere.center.y + sphere.radius * sinf(latNext),
          sphere.center.z + sphere.radius * cosf(latNext) * sinf(lon),
      };

      Vector3 c = {
          sphere.center.x + sphere.radius * cosf(lat) * cosf(lonNext),
          sphere.center.y + sphere.radius * sinf(lat),
          sphere.center.z + sphere.radius * cosf(lat) * sinf(lonNext),
      };

      // 視点変換＋NDC正規化
      Vector4 a4 = Transform4(a, viewProjectionMatrix);
      Vector4 b4 = Transform4(b, viewProjectionMatrix);
      Vector4 c4 = Transform4(c, viewProjectionMatrix);

      Vector3 ndcA = ToNDC(a4);
      Vector3 ndcB = ToNDC(b4);
      Vector3 ndcC = ToNDC(c4);

      // ビューポート変換
      Vector3 screenA = Transform(ndcA, viewportMatrix);
      Vector3 screenB = Transform(ndcB, viewportMatrix);
      Vector3 screenC = Transform(ndcC, viewportMatrix);

      Novice::DrawLine((int)screenA.x, (int)screenA.y, (int)screenB.x,
                       (int)screenB.y, color);
      Novice::DrawLine((int)screenA.x, (int)screenA.y, (int)screenC.x,
                       (int)screenC.y, color);
    }
  }
}

void DrawEdge(int index0, int index1, const Vector3 *screenVertices,
              uint32_t color) {
  Novice::DrawLine(static_cast<int>(screenVertices[index0].x),
                   static_cast<int>(screenVertices[index0].y),
                   static_cast<int>(screenVertices[index1].x),
                   static_cast<int>(screenVertices[index1].y), color);
}

void DrawAABB(const AABB &aabb, const Matrix4x4 &viewProjectionMatrix,
              const Matrix4x4 &viewportMatrix, uint32_t color) {

  // 8頂点を定義
  Vector3 vertices[8] = {
      {aabb.min.x, aabb.min.y, aabb.min.z},
      {aabb.max.x, aabb.min.y, aabb.min.z},
      {aabb.min.x, aabb.max.y, aabb.min.z},
      {aabb.max.x, aabb.max.y, aabb.min.z},
      {aabb.min.x, aabb.min.y, aabb.max.z},
      {aabb.max.x, aabb.min.y, aabb.max.z},
      {aabb.min.x, aabb.max.y, aabb.max.z},
      {aabb.max.x, aabb.max.y, aabb.max.z},
  };

  // スクリーン座標に変換
  Vector3 screenVertices[8];
  for (int i = 0; i < 8; ++i) {
    Vector4 clip;
    clip = Transform4(vertices[i], viewProjectionMatrix); // ViewProjection変換
    Vector3 ndc;
    ndc = ToNDC(clip); // NDCに正規化
    Vector3 screen;
    screen = Transform(ndc, viewportMatrix); // ビューポート変換
    screenVertices[i] = screen;
  }

  // 線を描画する関数

  // 底面
  DrawEdge(0, 1, screenVertices, color);
  DrawEdge(1, 3, screenVertices, color);
  DrawEdge(3, 2, screenVertices, color);
  DrawEdge(2, 0, screenVertices, color);

  // 天面
  DrawEdge(4, 5, screenVertices, color);
  DrawEdge(5, 7, screenVertices, color);
  DrawEdge(7, 6, screenVertices, color);
  DrawEdge(6, 4, screenVertices, color);

  // 側面
  DrawEdge(0, 4, screenVertices, color);
  DrawEdge(1, 5, screenVertices, color);
  DrawEdge(3, 7, screenVertices, color);
  DrawEdge(2, 6, screenVertices, color);
}

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

  // ライブラリの初期化
  Novice::Initialize(kWindowTitle, static_cast<int>(kWindowWidth),
                     static_cast<int>(kWindowHeight));

  // キー入力結果を受け取る箱
  char keys[256] = {0};
  char preKeys[256] = {0};

  Camera camera{
      {1.0f, 1.0f, 1.0f},  // scale
      {0.26f, 0.0f, 0.0f}, // rotate
      {0.0f, 1.9f, -6.49f} // translate
  };

  int wheel = 0;
  int mouseX = 0;
  int mouseY = 0;
  int prevMouseX = 0;
  int prevMouseY = 0;

  Vector3 scale{1.0f, 1.0f, 1.0f};
  Vector3 rotate{0.0f, 0.0f, 0.0f};
  Vector3 translate{0.0f, 0.0f, 0.0f};

  Sphere sphere = {{0.0f, 0.0f, 0.0f}, 1.0f};
  AABB aabb{.min{-0.5f, -0.5f, -0.5f}, .max{0.0f, 0.0f, 0.0f}};

  uint32_t color = WHITE;

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

    ImGui::DragFloat3("Sphere Center", &sphere.center.x, 0.01f);
    ImGui::DragFloat("Sphere Radius", &sphere.radius, 0.01f);
    ImGui::DragFloat3("AABB.min", &aabb.min.x, 0.01f);
    ImGui::DragFloat3("AABB.max", &aabb.max.x, 0.01f);

    ImGui::End();

    Novice::GetMousePosition(&mouseX, &mouseY);
    wheel = Novice::GetWheel();

    UpdateCameraControl(camera, keys, wheel, mouseX, mouseY, prevMouseX,
                        prevMouseY);

    Matrix4x4 worldMatrix = MakeAffineMatrix(scale, rotate, translate);
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

    //minとmaxが入れ替わらないようにする
    aabb.min.x = (std::min)(aabb.min.x, aabb.max.x);
    aabb.min.y = (std::min)(aabb.min.y, aabb.max.y);
    aabb.min.z = (std::min)(aabb.min.z, aabb.max.z);

    if (IsCollision(aabb, sphere)) {
      color = RED;
    } else {
      color = WHITE;
    }

    ///
    /// ↑更新処理ここまで
    ///

    ///
    /// ↓描画処理ここから
    ///

    DrawGrid(worldViewProjectionMatrix, viewportMatrix);

    DrawSphere(sphere, worldViewProjectionMatrix, viewportMatrix, color);

    DrawAABB(aabb, worldViewProjectionMatrix, viewportMatrix, color);

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