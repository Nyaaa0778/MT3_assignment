#include <Novice.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <imgui.h>
#include <stdint.h>

const char kWindowTitle[] = "LE2B_27_ヤマダ_ナオ_2_3_確認課題";

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

struct Segment {
  Vector3 origin; // 始点
  Vector3 diff;   // 終点への差分ベクトル
};

struct Plane {
  Vector3 normal;
  float distance;
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

// <summary>
/// 3次元ベクトルの正規化
/// </summary>
/// <param name="vector">正規化したいベクトル</param>
/// <returns>正規化されたベクトル</returns>
Vector3 Normalize(const Vector3 &vector) {
  float length =
      sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);

  if (length == 0.0f) {
    return {0.0f, 0.0f, 0.0f}; // 零ベクトルを返す
  }

  return {vector.x / length, vector.y / length, vector.z / length};
}

Vector3 Perpendicular(const Vector3 &v) {
  if (fabs(v.x) > fabs(v.z)) {
    return {-v.y, v.x, 0.0f};
  } else {
    return {0.0f, -v.z, v.y};
  }
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
/// 直線と平面の衝突判定
/// </summary>
/// <param name="segment"></param>
/// <param name="plane"></param>
/// <returns></returns>
bool IsCollision(const Segment &segment, const Plane &plane) {
  // 垂直判定
  float dot = Dot(plane.normal, segment.diff);

  // 垂直ならば平行なので衝突するわけがない
  if (dot == 0.0f) {
    return false;
  }

  // 媒介変数tを求める
  float t = (plane.distance - Dot(segment.origin, plane.normal)) /
            Dot(segment.diff, plane.normal);

  return (t >= 0.0f && t <= 1.0f);
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

/// <summary>
/// 平面の描画
/// </summary>
/// <param name="plane"></param>
/// <param name="viewProjectionMatrix"></param>
/// <param name="viewportMatrix"></param>
/// <param name="color"></param>
/// <param name="scale"></param>
void DrawPlane(const Plane &plane, const Matrix4x4 &viewProjectionMatrix,
               const Matrix4x4 &viewportMatrix, uint32_t color,
               float scale = 1.5f) {
  Vector3 center = {plane.normal.x * plane.distance,
                    plane.normal.y * plane.distance,
                    plane.normal.z * plane.distance};

  Vector3 perpendiculars[2];
  perpendiculars[0] = Normalize(Perpendicular(plane.normal));
  perpendiculars[1] = Normalize(Cross(plane.normal, perpendiculars[0]));

  // スケールを掛ける
  for (int i = 0; i < 2; ++i) {
    perpendiculars[i].x *= scale;
    perpendiculars[i].y *= scale;
    perpendiculars[i].z *= scale;
  }

  Vector3 corners[4];
  corners[0] = Add(center, Add(perpendiculars[0], perpendiculars[1]));
  corners[1] = Add(center, Add(perpendiculars[0],
                               {-perpendiculars[1].x, -perpendiculars[1].y,
                                -perpendiculars[1].z}));
  corners[2] = Add(
      center,
      Add({-perpendiculars[0].x, -perpendiculars[0].y, -perpendiculars[0].z},
          {-perpendiculars[1].x, -perpendiculars[1].y, -perpendiculars[1].z}));
  corners[3] = Add(center, Add({-perpendiculars[0].x, -perpendiculars[0].y,
                                -perpendiculars[0].z},
                               perpendiculars[1]));

  Vector3 screenPoints[4];
  for (int i = 0; i < 4; ++i) {
    Vector4 viewProj = Transform4(corners[i], viewProjectionMatrix);
    if (viewProj.w == 0.0f)
      return;
    Vector3 ndc = ToNDC(viewProj);
    screenPoints[i] = Transform(ndc, viewportMatrix);
  }

  for (int i = 0; i < 4; ++i) {
    int next = (i + 1) % 4;
    Novice::DrawLine((int)screenPoints[i].x, (int)screenPoints[i].y,
                     (int)screenPoints[next].x, (int)screenPoints[next].y,
                     color);
  }
}

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

  // ライブラリの初期化
  Novice::Initialize(kWindowTitle, 1280, 720);

  // キー入力結果を受け取る箱
  char keys[256] = {0};
  char preKeys[256] = {0};

  Vector3 cameraTranslate{0.0f, 1.9f, -6.49f};
  Vector3 cameraRotate{0.26f, 0.0f, 0.0f};

  Vector3 scale{1.0f, 1.0f, 1.0f};
  Vector3 rotate{0.0f, 0.0f, 0.0f};
  Vector3 translate{0.0f, 0.0f, 0.0f};

  Segment segment{{-2.0f, -1.0f, 0.0f}, {3.0f, 2.0f, 2.0f}};
  Plane plane = {{0.0f, 1.0f, 0.0f}, 0.0f};

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
    ImGui::DragFloat3("Camera Translate", &cameraTranslate.x, 0.01f);
    ImGui::DragFloat3("Camera Rotate", &cameraRotate.x, 0.01f);

    ImGui::DragFloat3("Plane.Normal", &plane.normal.x, 0.01f);
    ImGui::DragFloat("Plane.Distance", &plane.distance, 0.01f);
    plane.normal = Normalize(plane.normal);

    ImGui::DragFloat3("Segment.Origin", &segment.origin.x, 0.01f);
    ImGui::DragFloat3("Segment.Diff", &segment.diff.x, 0.01f);

    ImGui::End();

    if (keys[DIK_D]) {
      cameraTranslate.x += 0.03f;
    }

    if (keys[DIK_A]) {
      cameraTranslate.x -= 0.03f;
    }

    if (keys[DIK_W]) {
      cameraTranslate.z += 0.03f;
    }

    if (keys[DIK_S]) {
      cameraTranslate.z -= 0.03f;
    }

    Matrix4x4 worldMatrix = MakeAffineMatrix(scale, rotate, translate);
    Matrix4x4 cameraMatrix =
        MakeAffineMatrix({1.0f, 1.0f, 1.0f}, cameraRotate, cameraTranslate);
    Matrix4x4 viewMatrix = Inverse(cameraMatrix);
    Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(
        0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
    Matrix4x4 viewProjectionMatrix = Multiply(viewMatrix, projectionMatrix);
    Matrix4x4 worldViewProjectionMatrix =
        Multiply(worldMatrix, viewProjectionMatrix);
    Matrix4x4 viewportMatrix = MakeViewportMatrix(
        0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);

    Vector4 start = Transform4(segment.origin, viewProjectionMatrix);
    Vector4 end =
        Transform4(Add(segment.origin, segment.diff), viewProjectionMatrix);

    Vector3 ndcStart = ToNDC(start);
    Vector3 ndcEnd = ToNDC(end);

    Vector3 screenStart = Transform(ndcStart, viewportMatrix);
    Vector3 screenEnd = Transform(ndcEnd, viewportMatrix);

    IsCollision(segment, plane);

    ///
    /// ↑更新処理ここまで
    ///

    ///
    /// ↓描画処理ここから
    ///

    DrawPlane(plane, viewProjectionMatrix, viewportMatrix, 0x00FF00FF);

    DrawGrid(worldViewProjectionMatrix, viewportMatrix);

    if (IsCollision(segment, plane)) {
      Novice::DrawLine(
          static_cast<int>(screenStart.x), static_cast<int>(screenStart.y),
          static_cast<int>(screenEnd.x), static_cast<int>(screenEnd.y), RED);
    } else {
      Novice::DrawLine(
          static_cast<int>(screenStart.x), static_cast<int>(screenStart.y),
          static_cast<int>(screenEnd.x), static_cast<int>(screenEnd.y), WHITE);
    }

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