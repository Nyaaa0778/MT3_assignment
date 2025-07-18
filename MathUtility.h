#pragma once

struct Vector3 {
  float x;
  float y;
  float z;

  /// <summary>
  /// 和の複合演算子
  /// </summary>
  /// <param name="v"></param>
  /// <returns></returns>
  Vector3 &operator+=(const Vector3 &v) {
    x += v.x;
    y += v.y;
    z += v.z;

    return *this;
  }

  /// <summary>
  /// 差の複合演算子
  /// </summary>
  /// <param name="v"></param>
  /// <returns></returns>
  Vector3 &operator-=(const Vector3 &v) {

    x -= v.x;
    y -= v.y;
    z -= v.z;

    return *this;
  }

  /// <summary>
  /// 積の複合演算子
  /// </summary>
  /// <param name="v"></param>
  /// <returns></returns>
  Vector3 &operator*=(const Vector3 &v) {
    x *= v.x;
    y *= v.y;
    z *= v.z;

    return *this;
  }

  /// <summary>
  /// 商の複合演算子
  /// </summary>
  /// <param name="v"></param>
  /// <returns></returns>
  Vector3 &operator/=(const Vector3 &v) {
    x /= v.x;
    y /= v.y;
    z /= v.z;

    return *this;
  }
};

struct Matrix4x4 {
  float m[4][4];
};

/// <summary>
/// 行列の和
/// </summary>
/// <param name="m1"></param>
/// <param name="m2"></param>
/// <returns></returns>
Matrix4x4 Add(const Matrix4x4 &m1, const Matrix4x4 &m2);

/// <summary>
/// 行列の差
/// </summary>
/// <param name="m1"></param>
/// <param name="m2"></param>
/// <returns></returns>
Matrix4x4 Subtract(const Matrix4x4 &m1, const Matrix4x4 &m2);

/// <summary>
/// 行列の積
/// </summary>
/// <param name="m1"></param>
/// <param name="m2"></param>
/// <returns></returns>
Matrix4x4 Multiply(const Matrix4x4 &m1, const Matrix4x4 &m2);

/// <summary>
/// 3次元ベクトルの和
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
Vector3 Add(const Vector3 &v1, const Vector3 v2);

/// <summary>
/// 3次元ベクトルの差
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
Vector3 Subtract(const Vector3 v1, const Vector3 v2);

/// <summary>
/// 3次元ベクトルの積
/// </summary>
/// <param name="s"></param>
/// <param name="v"></param>
/// <returns></returns>
Vector3 Multiply(float s, const Vector3 &v);



/// ========================================
/// 演算子のオーバーロード
/// ========================================

/// <summary>
/// 3次元ベクトルの和
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
Vector3 operator+(const Vector3 &v1, const Vector3 &v2);

/// <summary>
/// 3次元ベクトルの差
/// </summary>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
Vector3 operator-(const Vector3 &v1, const Vector3 &v2);

/// <summary>
/// 3次元ベクトルの積
/// </summary>
/// <param name="v"></param>
/// <param name="s"></param>
/// <returns></returns>
Vector3 operator*(const Vector3 &v, float s);

/// <summary>
/// 3次元ベクトルの積
/// </summary>
/// <param name="v"></param>
/// <param name="s"></param>
/// <returns></returns>
Vector3 operator*(float s, const Vector3 &v);

/// <summary>
/// 3次元ベクトルの商
/// </summary>
/// <param name="v"></param>
/// <param name="s"></param>
/// <returns></returns>
Vector3 operator/(const Vector3 &v, float s);

/// <summary>
/// 行列の和
/// </summary>
/// <param name="m1"></param>
/// <param name="m2"></param>
/// <returns></returns>
Matrix4x4 operator+(const Matrix4x4 &m1, const Matrix4x4 &m2);

/// <summary>
/// 行列の差
/// </summary>
/// <param name="m1"></param>
/// <param name="m2"></param>
/// <returns></returns>
Matrix4x4 operator-(const Matrix4x4 &m1, const Matrix4x4 &m2);

/// <summary>
/// 行列の積
/// </summary>
/// <param name="m1"></param>
/// <param name="m2"></param>
/// <returns></returns>
Matrix4x4 operator*(const Matrix4x4 m1, const Matrix4x4 &m2);