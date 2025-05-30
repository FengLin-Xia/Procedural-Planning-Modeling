#define _MAIN_
#ifdef _MAIN_

#include "main.h"
#include <headers/zApp/include/zObjects.h>
#include <headers/zApp/include/zViewer.h>
#include <algorithm>
#include <unordered_map>

#include <array>
#include <vector>
#include <random>
#include <string>
#include <ctime>
#include <cmath>
#include <cfloat>
#include <set>
static std::mt19937 globalRng;
using namespace zSpace;


#ifndef M_PI
#define M_PI 3.14159265358979323846   // π 的双精度常量
#endif

/* ───────── CONSTANTS ───────── */
constexpr float mergeEPS = 5.0f;
constexpr float connectRadius = 20.0f;
constexpr float minEdge = 15.0f;
constexpr float angEPSdeg = 15.0f;

constexpr float sigma = 15.0f;   // 高斯核
constexpr int   gridN = 120;     // 向量场分辨率
zVector cpt{ -43.f, 12.f, 0.f };
/* ───────── HELPERS ───────── */
namespace zSpace {
    inline zVector operator+(const zVector& a, const zVector& b) { return zVector(a.x + b.x, a.y + b.y, a.z + b.z); }
    inline zVector operator-(const zVector& a, const zVector& b) { return zVector(a.x - b.x, a.y - b.y, a.z - b.z); }
    inline zVector operator*(const zVector& v, float s) { return zVector(v.x * s, v.y * s, v.z * s); }
}
inline Alice::vec z2A(const zVector& v) { return Alice::vec(v.x, v.y, v.z); }
inline float lerp(float a, float b, float t) { return a + (b - a) * t; }
inline float deg(float r) { return r * 57.2957795f; }

/* HSV→RGB */
static void hsv2rgb(float h, float s, float v, float& r, float& g, float& b)
{
    h = fmodf(h, 1.f) * 6.f; int i = int(h);
    float f = h - i, p = v * (1 - s), q = v * (1 - f * s), t = v * (1 - (1 - f) * s);
    switch (i) {
    case 0: r = v; g = t; b = p; break;
    case 1: r = q; g = v; b = p; break;
    case 2: r = p; g = v; b = t; break;
    case 3: r = p; g = q; b = v; break;
    case 4: r = t; g = p; b = v; break;
    default:r = v; g = p; b = q; break;
    }
}

/* ══════════════════════════════════════
   PART A  ——  样条网络 + 路网 / 向量场
   （与前版相同，这里直接给出实现）
   ══════════════════════════════════════ */

   /* ---------- SplineNet ---------- */
class SplineNet {
public:

    std::vector<std::pair<zVector, zVector>> rays;   // cpt → 最远三点的线
    std::vector<zVector> seed1;                      // 圆与直线交点
    std::vector<std::vector<zVector>> curves;        // 存圆

    void build(int samples = 60) {
        rays.clear();
        seed1.clear();
        curves.clear();

        // 1. 四个边界中点
        std::vector<zVector> borderCenters = {
            zVector(0, -50, cpt.z),  // bottom
            zVector(0, 50, cpt.z),   // top
            zVector(-50, 0, cpt.z),  // left
            zVector(50, 0, cpt.z)    // right
        };

        // 2. 排序并取最远的三个
        std::sort(borderCenters.begin(), borderCenters.end(), [&](const zVector& a, const zVector& b) {
            return lengthSquared(a - cpt) > lengthSquared(b - cpt);
            });

        for (int i = 0; i < 2; ++i) {
            zVector end = borderCenters[i];
            rays.emplace_back(cpt, end);

            zVector dir = normalize(end - cpt);
            zVector intersection = cpt + dir * 65.0f;
            seed1.push_back(intersection);

            std::cout
                << "seed" << i
                << " = ("
                << intersection.x << ", "
                << intersection.y << ", "
                << intersection.z << ")\n";
        }

        // 3. 绘制圆（半径30）
        std::vector<zVector> circle;
        for (int i = 0; i <= samples; ++i) {
            float t = float(i) / samples;
            float angle = t * 2.0f * M_PI;
            circle.emplace_back(
                cpt.x + 65.0f * std::cos(angle),
                cpt.y + 65.0f * std::sin(angle),
                cpt.z
            );
        }
        curves.push_back(std::move(circle));
    }

    void drawRays() const {
        glColor3f(0, 0, 0);
        for (const auto& r : rays) {
            drawLine(z2A(r.first), z2A(r.second));
        }
    }

    void drawCurves() const {
        glColor3f(0, 0, 0);
        if (!curves.empty()) {
            glBegin(GL_LINE_LOOP);
            for (const auto& pt : curves[0]) {
                glVertex3f(pt.x, pt.y, pt.z);
            }
            glEnd();
        }
    }

    void drawSeeds() const {
        glColor3f(1, 0, 0); // 红色点
        glPointSize(5.0f);
        glBegin(GL_POINTS);
        for (const auto& pt : seed1) {
            glVertex3f(pt.x, pt.y, pt.z);
        }
        glEnd();
    }


    // 在 SplineNet 类里添加这两个成员：
    std::vector<std::vector<zVector>> contour0Lines; // 存第1,3,5…层
    std::vector<std::vector<zVector>> contour1Lines; // 存第2,4,6…层
    std::vector<std::vector<zVector>> contour7Lines;  // 第七层的曲线集合


    // 然后用下面的定义替换原来的 drawContourMap：
    void drawContourMap(
        const zVector& seed0,
        const std::vector<zVector>& /*seed1*/,
        int gridRes /*=120*/,
        int nLevels /*=40*/,
        float spacing /*=0.25f*/
    ) {
        const float minX = -50.f, maxX = 50.f;
        const float minY = -50.f, maxY = 50.f;
        const float dx = (maxX - minX) / float(gridRes - 1);
        const float dy = (maxY - minY) / float(gridRes - 1);

        // 定义五个椭圆：包含 seed0 以及四个固定点，每个椭圆带 a、b、theta
        struct Ellipse { zVector c; float a, b, theta; };
        std::vector<Ellipse> ellipses;
        ellipses.reserve(5);
        // seed0 对应的椭圆
        float theta0 = 0.f;
        if (seed1.size() >= 2) {
            zVector mid = (seed1[0] + seed1[1]) * 0.5f;
            zVector vec = mid - seed0;
            theta0 = std::atan2(vec.y, vec.x);
        }
        ellipses.push_back({ seed0, 18.f, 15.f, theta0 });
        // 四个固定中心点
        for (int i = 0; i < seed1.size(); ++i) {
            const zVector& pt = seed1[i];
            const auto& ray = rays[i];

            zVector dir = ray.second - ray.first;
            float theta = std::atan2(dir.y, dir.x);  // 方向角（xy平面）

            // 添加椭圆，a=12长轴，b=6短轴，theta 为朝向方向
            ellipses.push_back({ pt, 10.f, 18.f, theta });
        }

        // 构建归一化半径场 D：对每个采样点，取到所有椭圆的最小归一化半径
        std::vector<std::vector<float>> D(gridRes, std::vector<float>(gridRes));
        auto evalNormRadius = [&](float x, float y) {
            float r_min = FLT_MAX;
            for (auto& e : ellipses) {
                // 旋转坐标，计算归一化半径
                float dx0 = x - e.c.x;
                float dy0 = y - e.c.y;
                float ux = dx0 * std::cos(-e.theta) - dy0 * std::sin(-e.theta);
                float uy = dx0 * std::sin(-e.theta) + dy0 * std::cos(-e.theta);
                float r_norm = std::sqrt((ux * ux) / (e.a * e.a) + (uy * uy) / (e.b * e.b));
                r_min = std::fmin(r_min, r_norm);
            }
            return r_min;
            };
        for (int i = 0; i < gridRes; ++i) {
            float x = minX + i * dx;
            for (int j = 0; j < gridRes; ++j) {
                float y = minY + j * dy;
                D[i][j] = evalNormRadius(x, y);
            }
        }

        // 清空历史数据
        contour0Lines.clear();
        contour1Lines.clear();

        // marching-squares 提取 iso = k*spacing 的归一化半径线
        glColor3f(0.4f, 0.4f, 0.4f);
        for (int k = 1; k <= nLevels; ++k) {
            float iso = spacing * k;
            std::vector<zVector> ptsBuf;
            for (int i = 0; i < gridRes - 1; ++i) {
                float x0 = minX + i * dx;
                float x1 = x0 + dx;
                for (int j = 0; j < gridRes - 1; ++j) {
                    float y0 = minY + j * dy;
                    float y1 = y0 + dy;
                    float v0 = D[i][j];
                    float v1 = D[i + 1][j];
                    float v2 = D[i + 1][j + 1];
                    float v3 = D[i][j + 1];
                    int m = (v0 > iso) | ((v1 > iso) << 1) | ((v2 > iso) << 2) | ((v3 > iso) << 3);
                    if (m == 0 || m == 15) continue;
                    auto interp = [&](float va, float vb, float a, float b) {
                        return a + (iso - va) / (vb - va) * (b - a);
                        };
                    zVector p[4];
                    if ((m & 1) != (m & 2)) p[0] = zVector(interp(v0, v1, x0, x1), y0, 0);
                    if ((m & 2) != (m & 4)) p[1] = zVector(x1, interp(v1, v2, y0, y1), 0);
                    if ((m & 8) != (m & 4)) p[2] = zVector(interp(v3, v2, x0, x1), y1, 0);
                    if ((m & 1) != (m & 8)) p[3] = zVector(x0, interp(v0, v3, y0, y1), 0);
                    auto drawSeg = [&](int a, int b) {
                        drawLine(z2A(p[a]), z2A(p[b]));
                        ptsBuf.push_back(p[a]); ptsBuf.push_back(p[b]);
                        };
                    switch (m) {
                    case 1: case 14: drawSeg(0, 3); break;
                    case 2: case 13: drawSeg(0, 1); break;
                    case 3: case 12: drawSeg(1, 3); break;
                    case 4: case 11: drawSeg(1, 2); break;
                    case 5: drawSeg(0, 1); drawSeg(3, 2); break;
                    case 6: case 9: drawSeg(0, 2); break;
                    case 7: case 8: drawSeg(2, 3); break;
                    }
                }
            }
            if (!ptsBuf.empty()) {
                if (k == 7) contour7Lines.emplace_back(ptsBuf);
                if (k % 2 == 1) contour0Lines.emplace_back(std::move(ptsBuf));
                else         contour1Lines.emplace_back(std::move(ptsBuf));
            }
        }

    }

    void drawExtendedCurves(
        const zVector& seed0,
        const std::vector<zVector>& seed1,
        const std::vector<std::pair<zVector, zVector>>& rays,
        std::vector<std::pair<zVector, zVector>>& yellowLines, // ✅ 收集线段
        int samples = 100)
    {
        const float minX = -50.f, maxX = 50.f;
        const float minY = -50.f, maxY = 50.f;

        auto inBox = [&](const zVector& pt) {
            return pt.x >= minX && pt.x <= maxX && pt.y >= minY && pt.y <= maxY;
            };

        auto extendTangentLine = [&](const zVector& pt, float angle) {
            zVector dir(std::cos(angle), std::sin(angle), 0);
            zVector A = pt;
            zVector B = pt + dir * 1000.f;

            zVector clipped;
            bool found = false;

            auto intersect = [&](float ex, float ey, bool vertical) {
                float u, ix, iy;
                if (vertical)
                    u = (ex - A.x) / (B.x - A.x), ix = ex, iy = A.y + (B.y - A.y) * u;
                else
                    u = (ey - A.y) / (B.y - A.y), iy = ey, ix = A.x + (B.x - A.x) * u;
                if (u >= 0.f && u <= 1.f && ix >= minX && ix <= maxX && iy >= minY && iy <= maxY) {
                    clipped = zVector(ix, iy, pt.z);
                    found = true;
                }
                };

            intersect(-50.f, 0.f, true);
            intersect(50.f, 0.f, true);
            intersect(0.f, -50.f, false);
            intersect(0.f, 50.f, false);

            if (found) {
                glColor3f(1.0f, 1.0f, 1.0f);
                drawLine(z2A(pt), z2A(clipped));
                yellowLines.emplace_back(pt, clipped);  // ✅ 收集切线
            }
            };

        auto drawHalfEllipse = [&](const zVector& center, float a, float b, float theta) -> std::vector<zVector> {
            std::vector<zVector> arc;
            for (int i = 0; i <= samples; ++i) {
                float t = float(i) / samples;
                float angle = -M_PI / 2.f + t * M_PI;
                float x = a * std::cos(angle);
                float y = b * std::sin(angle);
                float xr = x * std::cos(theta) - y * std::sin(theta);
                float yr = x * std::sin(theta) + y * std::cos(theta);
                arc.emplace_back(center.x + xr, center.y + yr, center.z);
            }

            glColor3f(1.0f, 1.0f, 1.0f);
            glBegin(GL_LINE_STRIP);
            for (auto& pt : arc)
                if (inBox(pt)) glVertex3f(pt.x, pt.y, pt.z);
            glEnd();

            for (size_t k = 0; k + 1 < arc.size(); ++k)
                yellowLines.emplace_back(arc[k], arc[k + 1]);  // ✅ 收集椭圆段

            return arc;
            };

        // 1. 大黄圆裁剪（保留画法）
        for (int i = 0; i < seed1.size(); ++i) {
            zVector dir = seed1[i] - rays[i].first;
            float len = std::sqrt(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);
            if (len < 1e-6f) continue;
            dir = dir * (1.0f / len);
            zVector center = seed1[i] + dir * 35.f;

            std::vector<zVector> segment;
            segment.reserve(samples + 1);
            bool inSeg = false;
            zVector prevPt;

            for (int s = 0; s <= samples; ++s) {
                float t = float(s) / samples;
                float a = t * 2.0f * M_PI;
                zVector pt{
                    center.x + 55.f * std::cos(a),
                    center.y + 55.f * std::sin(a),
                    center.z
                };

                bool inBox = (pt.x >= -50.f && pt.x <= 50.f &&
                    pt.y >= -50.f && pt.y <= 50.f);

                if (inBox) {
                    if (!inSeg && s > 0) {
                        float u, ix, iy;
                        auto testEdge = [&](float ex, float ey, bool vertical) {
                            if (vertical)
                                u = (ex - prevPt.x) / (pt.x - prevPt.x), ix = ex, iy = prevPt.y + (pt.y - prevPt.y) * u;
                            else
                                u = (ey - prevPt.y) / (pt.y - prevPt.y), iy = ey, ix = prevPt.x + (pt.x - prevPt.x) * u;
                            return u >= 0.f && u <= 1.f &&
                                ix >= -50.f && ix <= 50.f &&
                                iy >= -50.f && iy <= 50.f;
                            };
                        if (testEdge(-50.f, 0, true) ||
                            testEdge(50.f, 0, true) ||
                            testEdge(0, -50.f, false) ||
                            testEdge(0, 50.f, false))
                            segment.emplace_back(ix, iy, center.z);
                    }

                    segment.push_back(pt);
                    inSeg = true;
                }
                else if (inSeg) {
                    float u, ix, iy;
                    auto testEdge = [&](float ex, float ey, bool vertical) {
                        if (vertical)
                            u = (ex - prevPt.x) / (pt.x - prevPt.x), ix = ex, iy = prevPt.y + (pt.y - prevPt.y) * u;
                        else
                            u = (ey - prevPt.y) / (pt.y - prevPt.y), iy = ey, ix = prevPt.x + (pt.x - prevPt.x) * u;
                        return u >= 0.f && u <= 1.f &&
                            ix >= -50.f && ix <= 50.f &&
                            iy >= -50.f && iy <= 50.f;
                        };
                    if (testEdge(-50.f, 0, true) ||
                        testEdge(50.f, 0, true) ||
                        testEdge(0, -50.f, false) ||
                        testEdge(0, 50.f, false))
                        segment.emplace_back(ix, iy, center.z);

                    if (segment.size() >= 2) {
                        glColor3f(1.0f, 1.0f, 1.0f);
                        glBegin(GL_LINE_STRIP);
                        for (auto& p : segment)
                            glVertex3f(p.x, p.y, p.z);
                        glEnd();

                        for (size_t k = 0; k + 1 < segment.size(); ++k)
                            yellowLines.emplace_back(segment[k], segment[k + 1]);  // ✅ 收集圆段
                    }

                    segment.clear();
                    inSeg = false;
                }

                prevPt = pt;
            }

            if (inSeg && segment.size() >= 2) {
                glColor3f(1.0f, 1.0f, 1.0f);
                glBegin(GL_LINE_STRIP);
                for (auto& p : segment)
                    glVertex3f(p.x, p.y, p.z);
                glEnd();

                for (size_t k = 0; k + 1 < segment.size(); ++k)
                    yellowLines.emplace_back(segment[k], segment[k + 1]);  // ✅ 收集圆段
            }
        }

        // 2. seed0 半椭圆和两条切线
        float theta = 0.f;
        if (seed1.size() >= 2) {
            zVector mid = (seed1[0] + seed1[1]) * 0.5f;
            zVector vec = mid - seed0;
            theta = std::atan2(vec.y, vec.x);
        }

        std::vector<zVector> arc = drawHalfEllipse(seed0, 26.5f, 23.f, theta);
        if (arc.size() >= 2) {
            float angleLeft = theta + M_PI;
            float angleRight = theta + M_PI;
            extendTangentLine(arc.front(), angleLeft);
            extendTangentLine(arc.back(), angleRight);
        }

        // 3. seed0 → 最远边界延长线
        std::vector<zVector> borderCenters = {
            zVector(0, -50, seed0.z), zVector(0, 50, seed0.z),
            zVector(-50, 0, seed0.z), zVector(50, 0, seed0.z)
        };

        int imax = 0;
        float maxD2 = 0.0f;
        for (int i = 0; i < borderCenters.size(); ++i) {
            zVector d = borderCenters[i] - seed0;
            float d2 = d.x * d.x + d.y * d.y;
            if (d2 > maxD2) {
                maxD2 = d2;
                imax = i;
            }
        }

        zVector dir = borderCenters[imax] - seed0;
        float len = std::sqrt(dir.x * dir.x + dir.y * dir.y);
        if (len > 1e-6f) {
            dir = dir * (1.0f / len);
            zVector A = seed0 - dir * 1000.f;
            zVector B = seed0 + dir * 1000.f;

            std::vector<zVector> clipPts;
            auto intersect = [&](float ex, float ey, bool vertical) {
                float u, ix, iy;
                if (vertical)
                    u = (ex - A.x) / (B.x - A.x), ix = ex, iy = A.y + (B.y - A.y) * u;
                else
                    u = (ey - A.y) / (B.y - A.y), iy = ey, ix = A.x + (B.x - A.x) * u;
                if (u >= 0.f && u <= 1.f && ix >= minX && ix <= maxX && iy >= minY && iy <= maxY)
                    clipPts.emplace_back(ix, iy, seed0.z);
                };

            intersect(-50.f, 0.f, true);
            intersect(50.f, 0.f, true);
            intersect(0.f, -50.f, false);
            intersect(0.f, 50.f, false);

            if (clipPts.size() == 2) {
                glColor3f(1.0f, 1.0f, 1.0f);
                drawLine(z2A(clipPts[0]), z2A(clipPts[1]));
                yellowLines.emplace_back(clipPts[0], clipPts[1]);  // ✅ 收集主线
            }
        }
    }






    std::vector<zVector> computeCircleIntersections() const {
        return seed1;

    }

private:
    static float lengthSquared(const zVector& v) {
        return v.x * v.x + v.y * v.y + v.z * v.z;
    }

    static zVector normalize(const zVector& v) {
        float len = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
        if (len < 1e-6f) return zVector(0, 0, 0);
        return zVector(v.x / len, v.y / len, v.z / len);
    }
    static float dot(const zVector& a, const zVector& b) {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

};



/* ---------- MarkSet : 白 / 黄 / 青点 ---------- */
struct MarkSet {
    std::vector<zVector> merged;
    std::vector<zVector> intPts;

    static bool segIntersect2D(const zVector& a, const zVector& b,
        const zVector& c, const zVector& d,
        zVector& out)
    {
        float x1 = a.x, y1 = a.y, x2 = b.x, y2 = b.y;
        float x3 = c.x, y3 = c.y, x4 = d.x, y4 = d.y;
        float den = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1);
        if (fabs(den) < 1e-6f) return false;
        float ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / den;
        float ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / den;
        if (ua < 0 || ua > 1 || ub < 0 || ub > 1) return false;
        out = zVector(x1 + ua * (x2 - x1), y1 + ua * (y2 - y1), 0);
        return true;
    }

    void build(const SplineNet& net)
    {
        // 只收集交点
        struct Raw { zVector p; };
        std::vector<Raw> raw;

        // 射线-射线交点
        for (size_t i = 0; i < net.rays.size(); ++i) {
            for (size_t j = i + 1; j < net.rays.size(); ++j) {
                zVector P;
                if (segIntersect2D(net.rays[i].first, net.rays[i].second,
                    net.rays[j].first, net.rays[j].second, P))
                {
                    raw.push_back({ P });
                }
            }
        }
        // 射线-样条交点
        for (auto& ray : net.rays) {
            for (auto& C : net.curves) {
                for (size_t k = 0; k + 1 < C.size(); ++k) {
                    zVector P;
                    if (segIntersect2D(ray.first, ray.second, C[k], C[k + 1], P)) {
                        raw.push_back({ P });
                    }
                }
            }
        }



        // 简单聚类：将靠近的交点合并
        struct Cluster { zVector ip; };
        std::vector<Cluster> cls;
        for (auto& m : raw) {
            bool fit = false;
            for (auto& c : cls) {
                if ((m.p - c.ip).length() < mergeEPS) {
                    // 平均到已有中心
                    c.ip = (c.ip + m.p) * 0.5f;
                    fit = true;
                    break;
                }
            }
            if (!fit) {
                cls.push_back({ m.p });
            }
        }

        // 输出
       // 输出，跳过 x==50 或 y==50 的点
        merged.clear();
        intPts.clear();
        for (auto& c : cls) {
            // 如果横坐标或纵坐标正好是 50，就不加入 seed1
            if (c.ip.x == 50.0f || c.ip.y == 50.0f) continue;

            merged.push_back(c.ip);
            intPts.push_back(c.ip);
        }

    }




    void draw() const {
        // 只画交点簇中心
        glColor3f(0, 0, 0);
        for (auto& p : intPts)
            drawCircle(z2A(p), 1.4, 16);
    }



};





/* ---------- EdgeSet : 洋红连线 + 橙中点 ---------- */





/* ---------- VectorField (对比度更大) ---------- */
class VectorField {
    float cell, minV, maxV;
    float scalar[gridN][gridN];
public:
    void compute(const std::vector<zVector>& src, const std::vector<float>& w)
    {
        float half = 50.f; cell = (half * 2) / (gridN - 1);
        minV = 1e9f; maxV = -1e9f;
        for (int i = 0; i < gridN; ++i)for (int j = 0; j < gridN; ++j) {
            float x = -half + i * cell, y = -half + j * cell, v = 0;
            for (size_t k = 0; k < src.size(); ++k) {
                float d = sqrtf((x - src[k].x) * (x - src[k].x) + (y - src[k].y) * (y - src[k].y));
                v += w[k] * 0.3f * expf(-(d * d) / (sigma * sigma));      // ×0.3
            }
            scalar[i][j] = v; minV = fminf(minV, v); maxV = fmaxf(maxV, v);
        }
    }
    void draw()const {
        float half = 50.f;
        for (int i = 0; i < gridN - 1; ++i)for (int j = 0; j < gridN - 1; ++j) {
            float v[4] = { scalar[i][j],scalar[i + 1][j],
                        scalar[i + 1][j + 1],scalar[i][j + 1] };
            float xs[4] = { -half + i * cell,-half + (i + 1) * cell,
                         -half + (i + 1) * cell,-half + i * cell };
            float ys[4] = { -half + j * cell,-half + j * cell,
                         -half + (j + 1) * cell,-half + (j + 1) * cell };
            glBegin(GL_QUADS);
            for (int k = 0; k < 4; ++k) {
                float t = (v[k] - minV) / (maxV - minV + 1e-6f);
                t = powf(t, 0.6f);                             // γ-校正
                float r, g, b; hsv2rgb(0.55f + 0.4f * t, 1, 1, r, g, b);
                glColor3f(r, g, b); glVertex3f(xs[k], ys[k], 0);
            }
            glEnd();
        }
    }






};

/* ══════════════════════════════════════
   PART B  ——  ColorGrowth class
   ══════════════════════════════════════ */

#include "ColorGrowth.h"  // 定义 ColorType、GPoint



   // ---------- EllipseExpansion.h ----------

#include "EllipseExpansion.h"



/* ══════════════════════════════════════
   PART C  ——  Scene : 5-stage state machine
   ══════════════════════════════════════ */

struct Snapshot {
    std::vector<std::vector<zVector>>       curves;
    std::vector<std::pair<zVector, zVector>> rays;
    std::vector<GPoint>                     points;
    std::vector<GLine>                      lines;
};

class Scene {
    // —— 新增成员 —— 
    uint32_t      currentSeed;   // 本次流程用的种子
    std::string   cmdLog;        // 按键历史
    std::unordered_map<uint32_t, Snapshot> snapshots;
    std::vector<std::pair<zVector, zVector>> yellowLines;
    // —— 原有成员 —— 
    int           stage = 0;
    bool          showField = false;
    bool          showEllipseMode = false;
    bool showTriangleOnly = false;
    SplineNet        net;
    MarkSet          marks;

    VectorField      vField;
    ColorGrowth      growth;
    EllipseExpansion ellipses;
    bool showContours = false;

    bool    countRequested = false;
    size_t  prevGrowthCount = 0;

public:
    static constexpr int SMOOTH_PASSES = 2;

    Scene()
        : currentSeed(static_cast<uint32_t>(std::time(nullptr)))
    {
        globalRng.seed(currentSeed);
        cmdLog.clear();

        rebuild();
    }

    void Scene::rebuild() {
        globalRng.seed(currentSeed);
        showEllipseMode = false;

        net.drawExtendedCurves(::cpt, net.seed1, net.rays, yellowLines);

        // 1) 重建网格、交点、边网络
        net.build();
        marks.build(net);

        yellowLines.clear();
        // 2) 重置 growth，并添加初始种子
        growth.clear();
        growth.setPlacementMode(ColorGrowth::PlacementMode::Contour
        );
        growth.addSeed(0, cpt);

        // 3) （可选）把圆弧交点也塞进 ColorGrowth，如果你还要保留这部分逻辑
        {
            std::vector<zVector> circPts = net.computeCircleIntersections();
            std::vector<Alice::vec> apts;
            for (auto& z : circPts) apts.emplace_back(z.x, z.y, z.z);
            growth.setCircleIntersectionPoints(apts);
        }

        // 4) 清空并生成第1、3、5…条等高线到 net.contour0Lines
        net.contour0Lines.clear();
        net.contour1Lines.clear();
        net.drawContourMap(cpt, marks.intPts, 120, 40, 0.25);

        // 5) 把这些“0类”等高线线段传给 growth
        growth.setContour0Lines(net.contour0Lines);
        ellipses.setContour1Lines(net.contour1Lines);

        // 6) 重置状态机阶段
        stage = 0;
    }



    void next() {
        if (showEllipseMode) return;
        stage = (stage + 1) % 3;
        if (stage == 1)    for (auto& p : marks.merged) growth.addSeed(1, p);

    }

    void draw() {

        // —— T 模式优先：只画三角形 —— 
        if (showTriangleOnly) {
            ellipses.drawVectorTriangles();
            return;
        }

        if (showEllipseMode) {
            // 1) 先画等高线
           // net.drawContourMap(cpt, marks.intPts, 120, 40, 0.25f);
            // 2) 更新并画椭圆
            ellipses.update();
            ellipses.draw();

            // 3) 最后画黄线
            glColor3f(1.0f, 1.0f, 1.0f);  // 黄色
            for (auto& seg : yellowLines) {
                drawLine(z2A(seg.first), z2A(seg.second));
            }

            return;
        }
        if (showField) {
            std::vector<zVector> src; std::vector<float> w;
            src.push_back(cpt); w.push_back(0.3f);
            if (stage >= 1) for (auto& p : marks.merged) src.push_back(p), w.push_back(0.18f);

            vField.compute(src, w);
            vField.draw();
            return;
        }




        net.drawContourMap(cpt, marks.intPts, 120, 40, 0.25);
        net.drawExtendedCurves(::cpt, net.seed1, net.rays, yellowLines);
        // net.drawCurves();
        // net.drawRays();
        marks.draw();

        size_t before = prevGrowthCount;
        growth.update();
        size_t after = growth.getPoints().size();
        if (countRequested && before == after) {
            std::cout << "Total points after full-expand: " << after << std::endl;
            countRequested = false;
        }
        prevGrowthCount = after;
        growth.draw();

        if (stage == 2) {

            std::vector<Alice::vec> ctr;

            growth.recolorNear(ctr, 1.0f);
        }


    }

    // 保存当前完整状态到 snapshots[currentSeed]
    void recordSnapshot() {
        Snapshot s;
        s.curves = net.curves;
        s.rays = net.rays;
        s.points = growth.getPoints();
        s.lines = growth.getLines();
        snapshots[currentSeed] = std::move(s);
    }

    // 从快照恢复状态
    void loadSnapshot(uint32_t seed) {
        auto it = snapshots.find(seed);
        if (it == snapshots.end()) {
            std::cerr << "No snapshot for seed " << seed << std::endl;
            return;
        }
        const Snapshot& s = it->second;
        net.curves = s.curves;
        net.rays = s.rays;
        growth.loadState(s.points, s.lines);
        showEllipseMode = false;
        stage = 2;  // 恢复到“生长完毕”阶段
    }

    void handleKey(unsigned char k, int x, int y) {
        switch (k) {
        case 'v': case 'V':
            showField = !showField;
            break;
        case 'r': case 'R':
            currentSeed = static_cast<uint32_t>(std::time(nullptr));
            rebuild();
            cmdLog += 'r';
            break;
        case '1':
            growth.expandTo50();
            cmdLog += '1';
            break;
        case '2':
            growth.expandFull();
            countRequested = true;
            prevGrowthCount = growth.getPoints().size();
            cmdLog += '2';
            break;
        case 'a': case 'A':
            if (showEllipseMode) {
                ellipses.cornerSmooth(SMOOTH_PASSES);
                cmdLog += 'a';
            }
            break;
        case 't': case 'T':

            // 切换“仅三角形”开关
            showTriangleOnly = !showTriangleOnly;
            // 通知 EllipseExpansion 切换内部的 vector+triangle 开关
            ellipses.toggleVectorTriangle();
            cmdLog += 't';
            // 请求重绘
            glutPostRedisplay();

            break;
        case 'z': {
            showEllipseMode = true;
            cmdLog += 'z';

            // 先把上一轮的黄线清空
            yellowLines.clear();
            // 只在进入椭圆模式时收集一次
            net.drawExtendedCurves(cpt, net.seed1, net.rays, yellowLines);

            // 1) 构造样条线段列表
            std::vector<std::pair<zVector, zVector>> sl;
            for (auto& C : net.curves) {
                for (size_t i = 0; i + 1 < C.size(); ++i) {
                    sl.emplace_back(C[i], C[i + 1]);
                }
            }
            // 2) 拿到当前所有 growth 种子
            const auto gpts = growth.getPoints();  // vector<GPoint>
            // 3) edgeLines（如果你之前有收集的话），否则就传空
            std::vector<std::pair<zVector, zVector>> el;

            // 4) 初始化椭圆扩散，把收集到的黄线一次性传进去
            ellipses.init(
                cpt,               // seed0
                marks.intPts,      // seed1
                gpts,              // growthPts
                sl,                // spline 线段
                net.rays,          // ray 线段
                el,                // edge 线段
                yellowLines        // 刚刚收集的黄线
            );

            // 5) 设置等高线
            net.drawContourMap(cpt, marks.intPts, 120, 40, 0.25f);
            ellipses.setContour1Lines(net.contour1Lines);
            break;
        }

        case 's': case 'S':
            if (showEllipseMode) {
                ellipses.startDeform();
                cmdLog += 's';
            }
            break;
        case 'c': case 'C': {
            cmdLog += 'c';

            // 把 yellowLines（zVector 对）转成 GLine 数组
            std::vector<GLine> yln;
            yln.reserve(yellowLines.size());
            for (auto& seg : yellowLines) {
                yln.emplace_back(z2A(seg.first), z2A(seg.second));
            }


            growth.removeNearNetwork(yln, 1.5f);

            recordSnapshot();
            break;
        }
        case 'm': case 'M':
            // 删除所有距紫色点半径 12 单位内的异色点
            growth.removePointsNearPurple(6.0f);
            cmdLog += 'm';
            break;
        case 'n': case 'N':
            // 删除所有距紫色点半径 12 单位内的异色点
            growth.removePointsNearYellow(4.0f);
            cmdLog += 'n';
            break;
        case 'o': case 'O':
            if (showEllipseMode) {
                ellipses.toggleOffset();
                cmdLog += 'o';
            }
            break;
        case '0':
            std::cout << "ReplaySeed=" << currentSeed << " cmds=" << cmdLog << std::endl;
            break;
        case 'l': case 'L': {
            uint32_t seed;
            std::cout << "请输入 Seed 并回车：";
            if (std::cin >> seed) {
                loadSnapshot(seed);
            }
            else {
                std::cin.clear();
                std::string dummy; std::getline(std::cin, dummy);
                std::cerr << "Seed 输入失败，取消恢复\n";
            }
            break;
        }
        default:
            break;
        }
    }
};


/* ───────── GLOBAL Scene ───────── */
static Scene scene;

/* ───────── APP HOOKS ───────── */
void setup() {}
void update(int) {}
void draw() { backGround(1.0f);  scene.draw(); }
void mousePress(int b, int st, int, int) { if (b == 0 && st == 0) scene.next(); }
void keyPress(unsigned char k, int x, int y) {
    scene.handleKey(k, x, y);
}
void mouseMotion(int, int) {}

#endif // _MAIN_
