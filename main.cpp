#include <Novice.h>
#define _USE_MATH_DEFINES
#include <imgui.h>
#include <stdint.h>
#include"AffineMatrix.h"

const char kWindowTitle[] = "LE2B_27_ヤマダ_ナオ_3_0_確認課題";

const float kWindowWidth = 1280.0f;
const float kWindowHeight = 720.0f;

struct Vector4 {
  float x;
  float y;
  float z;
  float w;
};

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

/// <summary>
/// 距離の計算
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
float Length(const Vector3 &v1, const Vector3 &v2) {
  return sqrtf(powf(v1.x - v2.x, 2) + powf(v1.y - v2.y, 2) +
               powf(v1.z - v2.z, 2));
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
/// 線形補間
/// </summary>
/// <param name="p0"></param>
/// <param name="p1"></param>
/// <param name="t"></param>
/// <returns></returns>
Vector3 Lerp(const Vector3 &p0, const Vector3 &p1, float t) {
  Vector3 result;
  result.x = p0.x * (1.0f - t) + p1.x * t;
  result.y = p0.y * (1.0f - t) + p1.y * t;
  result.z = p0.z * (1.0f - t) + p1.z * t;

  return result;
}

/// <summary>
/// ベジエ曲線の式
/// </summary>
/// <param name="p0"></param>
/// <param name="p1"></param>
/// <param name="p2"></param>
/// <param name="t"></param>
/// <returns></returns>
Vector3 Bezier(const Vector3 &p0, const Vector3 &p1, const Vector3 &p2,
               float t) {
  Vector3 p0p1 = Lerp(p0, p1, t);
  Vector3 p1p2 = Lerp(p1, p2, t);

  return Lerp(p0p1, p1p2, t);
}

/// <summary>
/// ベジエ曲線を描画
/// </summary>
/// <param name="p0"></param>
/// <param name="p1"></param>
/// <param name="p2"></param>
/// <param name="worldViewProj"></param>
/// <param name="viewport"></param>
/// <param name="color"></param>
void DrawBezier(const Vector3 &p0, const Vector3 &p1, const Vector3 &p2,
                const Matrix4x4 &worldViewProj, const Matrix4x4 &viewport,
                uint32_t color) {

  // t=0 の点（始点）をサンプリング
  Vector3 prev = p0;
  Vector4 prev4 = Transform4(prev, worldViewProj);
  Vector3 prevN = ToNDC(prev4);
  Vector3 prevScreen = Transform(prevN, viewport);

  // 分割して線を描く
  for (int i = 1; i <= 32; ++i) {
    float t = i / float(32);
    Vector3 cur = Bezier(p0, p1, p2, t);

    // カレント位置をスクリーンに変換
    Vector4 cur4 = Transform4(cur, worldViewProj);
    if (cur4.w == 0.0f)
      continue; // w=0 だと正規化できないのでスキップ

    Vector3 curN = ToNDC(cur4);
    Vector3 curScreen = Transform(curN, viewport);

    // 線を描く
    Novice::DrawLine(int(prevScreen.x), int(prevScreen.y), int(curScreen.x),
                     int(curScreen.y), color);

    prevScreen = curScreen;
  }

  // 始点 p0 に黒丸を描く
  Vector4 p0_4 = Transform4(p0, worldViewProj);
  if (p0_4.w != 0.0f) {
    Vector3 p0_ndc = ToNDC(p0_4);
    Vector3 p0_screen = Transform(p0_ndc, viewport);
    Novice::DrawEllipse(int(p0_screen.x), int(p0_screen.y), 4, 4, 0.0f, BLACK,
                        kFillModeSolid);
  }

  // 終点 p2 に黒丸を描く
  Vector4 p2_4 = Transform4(p2, worldViewProj);
  if (p2_4.w != 0.0f) {
    Vector3 p2_ndc = ToNDC(p2_4);
    Vector3 p2_screen = Transform(p2_ndc, viewport);
    Novice::DrawEllipse(int(p2_screen.x), int(p2_screen.y), 4, 4, 0.0f, BLACK,
                        kFillModeSolid);
  }
}

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

  // ライブラリの初期化
  Novice::Initialize(kWindowTitle, static_cast<int>(kWindowWidth),
                     static_cast<int>(kWindowHeight));

  // キー入力結果を受け取る箱
  char keys[256] = {0};
  char preKeys[256] = {0};

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

  Vector3 controlPoints[3] = {
      {-0.8f, 0.58f, 1.0f}, {1.76f, 1.0f, -0.3f}, {0.94f, -0.7f, 2.3f}};

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

    ImGui::Text("ControlPoints");
    ImGui::DragFloat3("ControlPoint 0", &controlPoints[0].x, 0.01f);
    ImGui::DragFloat3("ControlPoint 1", &controlPoints[1].x, 0.01f);
    ImGui::DragFloat3("ControlPoint 2", &controlPoints[2].x, 0.01f);

    ImGui::End();

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

    DrawBezier(controlPoints[0], controlPoints[1], controlPoints[2],
               worldViewProjectionMatrix, viewportMatrix, BLUE);

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