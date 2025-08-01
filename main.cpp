#include "../../../KamataEngine/DirectXGame/math/Vector4.h"
#include "AffineMatrix.h"
#include <Novice.h>
#include <algorithm>
#include <imgui.h>

const char kWindowTitle[] = "LE2B_27_ヤマダ_ナオ_タ4_4_確認課題";

struct Plane {
  Vector3 normal;
  float distance;
};

struct Ball {
  Vector3 pos;
  Vector3 velocity;
  Vector3 acceleration;
  float mass;
  float radius;
  unsigned int color;
};

struct Sphere {
  Vector3 center;
  float radius;
};

struct Transformation {
  Vector3 scale;
  Vector3 rotate;
  Vector3 translate;
};

struct Segment {
  Vector3 origin;   // 始点
  Vector3 terminus; // 終点
};

struct Capsule {
  Segment segment; // 中心軸の線分
  float radius;    // 両端の球と側面の半径
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
  if (keys[DIK_D]) {
    camera.translate.x += 0.03f;
  }

  if (keys[DIK_A]) {
    camera.translate.x -= 0.03f;
  }

  if (keys[DIK_W]) {
    camera.translate.z += 0.03f;
  }

  if (keys[DIK_S]) {
    camera.translate.z -= 0.03f;
  }

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

/// <summary>
/// 反射ベクトルを求める
/// </summary>
/// <param name="input"></param>
/// <param name="normal"></param>
/// <returns></returns>
Vector3 Reflect(const Vector3 &input, const Vector3 &normal) {
  return input - 2.0f * Dot(input, normal) * normal;
}

/// <summary>
/// 正射影ベクトルを求める
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
Vector3 Project(const Vector3 &v1, const Vector3 &v2) {
  return Dot(v1, v2) /
         (std::powf(v2.x, 2) + std::powf(v2.y, 2) + std::powf(v2.z, 2)) * v2;
}

/// <summary>
/// ボールと平面の当たり判定
/// </summary>
/// <param name="sphere"></param>
/// <param name="plane"></param>
/// <returns></returns>
bool IsCollision(const Sphere &sphere, const Plane &plane) {
  // 平面の法線を正規化
  Vector3 n = Normalize(plane.normal);

  // 球の中心と平面との符号付き距離を計算
  //    plane.distance は「dot(n, 平面上の任意の点)」で事前に求めておく
  float signedDist = Dot(n, sphere.center) - plane.distance;

  // 絶対値が半径以下なら衝突
  return fabsf(signedDist) <= sphere.radius;
}

Vector3 Perpendicular(const Vector3 &v) {
  if (fabs(v.x) > fabs(v.z)) {
    return {-v.y, v.x, 0.0f};
  } else {
    return {0.0f, -v.z, v.y};
  }
}

bool IsCollision(const Capsule &capsule, const Plane &plane) {
  Vector3 n = Normalize(plane.normal);

  float d0 = fabsf(Dot(n, capsule.segment.origin) - plane.distance);
  float d1 = fabsf(Dot(n, capsule.segment.terminus) - plane.distance);

  // min(d0, d1) の代わりに if 文で判定
  float closest = (d0 < d1) ? d0 : d1;

  return closest <= capsule.radius;
}

/// <summary>
/// 平面の描画
/// </summary>
/// <param name="plane"></param>
/// <param name="viewProjectionMatrix"></param>
/// <param name="viewportMatrix"></param>
/// <param name="color"></param>
void DrawPlane(const Plane &plane, const Matrix4x4 &viewProjectionMatrix,
               const Matrix4x4 &viewportMatrix, uint32_t color = WHITE) {
  Vector3 center = Multiply(plane.distance, plane.normal);

  // 正しい2つの直交ベクトルを作る
  Vector3 tangent = Normalize(Perpendicular(plane.normal));
  Vector3 bitangent = Normalize(Cross(plane.normal, tangent));

  // 面積サイズ（2.0f で伸ばしている）
  float halfSize = 2.0f;

  Vector3 corners[4];
  corners[0] = Add(
      center, Add(Multiply(halfSize, tangent), Multiply(halfSize, bitangent)));
  corners[1] = Add(
      center, Add(Multiply(-halfSize, tangent), Multiply(halfSize, bitangent)));
  corners[2] = Add(center, Add(Multiply(-halfSize, tangent),
                               Multiply(-halfSize, bitangent)));
  corners[3] = Add(
      center, Add(Multiply(halfSize, tangent), Multiply(-halfSize, bitangent)));

  for (int i = 0; i < 4; ++i) {
    Vector4 clip = Transform4(corners[i], viewProjectionMatrix);
    if (clip.w != 0.0f) {
      Vector3 ndc = ToNDC(clip);
      corners[i] = Transform(ndc, viewportMatrix);
    }
  }

  for (int i = 0; i < 4; ++i) {
    int next = (i + 1) % 4;
    Novice::DrawLine(static_cast<int>(corners[i].x),
                     static_cast<int>(corners[i].y),
                     static_cast<int>(corners[next].x),
                     static_cast<int>(corners[next].y), color);
  }
}

/// <summary>
/// ボールの描画
/// </summary>
/// <param name="sphere"></param>
/// <param name="viewProjectionMatrix"></param>
/// <param name="viewportMatrix"></param>
/// <param name="color"></param>
void DrawSphere(const Sphere &sphere, const Matrix4x4 &viewProjectionMatrix,
                const Matrix4x4 &viewportMatrix, uint32_t color = WHITE) {
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

static int kWindowWidth = 1280;
static int kWindowHeight = 720;

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

  // ライブラリの初期化
  Novice::Initialize(kWindowTitle, 1280, 720);

  // キー入力結果を受け取る箱
  char keys[256] = {0};
  char preKeys[256] = {0};

  Plane plane = {
      Normalize({-0.2f, 0.9f, -0.3f}), // 法線
      0.0f                             // 距離
  };

  Ball ball = {
      {0.8f, 1.2f, 0.3f},  // 位置
      {0.0f, 0.0f, 0.0f},  // 速度
      {0.0f, -9.8f, 0.0f}, // 加速度
      2.0f,                // 質量
      0.05f,               // 半径
      WHITE                // 色
  };

  float deltaTime = 1.0f / 60.0f;
  float e = 0.5f;

  int isStart = false;

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
      isStart = true;
    }

    ImGui::End();

    Vector3 prevBallPos = ball.pos; // 追加：前フレームの位置

    if (isStart) {
      prevBallPos = ball.pos; // 前回位置を記録

      // 速度更新
      ball.velocity += ball.acceleration * deltaTime;

      // 位置更新
      ball.pos += ball.velocity * deltaTime;

      // 通常の球と平面の当たり判定
      if (IsCollision(Sphere{ball.pos, ball.radius}, plane)) {
        Vector3 reflected = Reflect(ball.velocity, plane.normal);
        Vector3 projectToNormal = Project(reflected, plane.normal);
        Vector3 movingDirection = reflected - projectToNormal;
        ball.velocity = projectToNormal * e + movingDirection;
      } else {
        // --- カプセルすり抜けチェック ---
        Capsule capsule = {{prevBallPos, ball.pos}, ball.radius};

        if (IsCollision(capsule, plane)) {
          // カプセルが接触してる = すり抜けた
          // → 平面に戻す補正
          Vector3 center = (prevBallPos + ball.pos) * 0.5f;
          float dist = Dot(plane.normal, center) - plane.distance;

          // 平面の外へ押し戻す（仮補正）
          ball.pos -= plane.normal * (dist - ball.radius);

          // 反射処理
          Vector3 reflected = Reflect(ball.velocity, plane.normal);
          Vector3 projectToNormal = Project(reflected, plane.normal);
          Vector3 movingDirection = reflected - projectToNormal;
          ball.velocity = projectToNormal * e + movingDirection;
        }
      }
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
    DrawPlane(plane, worldViewProjectionMatrix, viewportMatrix);
    DrawSphere(Sphere{ball.pos, ball.radius}, worldViewProjectionMatrix,
               viewportMatrix);

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
