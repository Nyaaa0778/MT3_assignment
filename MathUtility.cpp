#include "MathUtility.h"
#include<cmath>

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
/// 3次元ベクトルの和
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
/// 3次元ベクトルの差
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
Vector3 Subtract(const Vector3 v1, const Vector3 &v2) {
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

// <summary>
/// 3次元ベクトルの正規化
/// </summary>
/// <param name="vector">正規化したいベクトル</param>
/// <returns>正規化されたベクトル</returns>
Vector3 Normalize(const Vector3 &v) {
  float length = Length(v);

  if (length == 0.0f) {
    return {0.0f, 0.0f, 0.0f}; // 零ベクトルを返す
  }

  return v / length;
}

/// <summary>
/// 長さ
/// </summary>
/// <param name="v"></param>
/// <returns></returns>
float Length(const Vector3 &v) {
  return std::sqrtf(std::powf(v.x, 2) + std::powf(v.y, 2) + std::powf(v.z, 2));
}


/// <summary>
/// 内積
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
float Dot(const Vector3& v1, const Vector3& v2) {
  return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

/// <summary>
/// 外積
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
Vector3 Cross(const Vector3& v1, const Vector3& v2) {
  return {v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
          v1.x * v2.y - v1.y * v2.x};
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