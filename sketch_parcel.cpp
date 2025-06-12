#define _MAIN_
#ifdef _MAIN_

#include "main.h"
#include <headers/zApp/include/zObjects.h>
#include <headers/zApp/include/zViewer.h>
#include <algorithm>
#include <unordered_map>
#include <fstream>    // std::ifstream
#include <sstream>    // std::stringstream
#include <iostream>   // std::cerr, std::cout
#include <utility>    // std::pair
#include <filesystem> // C++17 的 std::filesystem::current_path()
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
zVector cpt{ -48.f, -18.f, 0.f };
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


std::vector<zVector> loadTerrainPoints(const std::string& filename, int stepX = 3, int stepY = 3)
{
    std::vector<zVector> pts;
    std::ifstream in(filename);
    if (!in.is_open()) {
        std::cerr << "Cannot open " << filename << "\n";
        return pts;
    }

    std::string line;
    std::vector<std::vector<zVector>> rows; // 按y分组，每一排
    float lastY = std::numeric_limits<float>::max();
    std::vector<zVector> currentRow;

    // 先逐行读取，按y分行
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string xstr, ystr;
        if (std::getline(ss, xstr, ',') && std::getline(ss, ystr, ',')) {
            float x = std::stof(xstr);
            float y = std::stof(ystr);

            if (currentRow.empty() || std::fabs(y - lastY) < 1e-4) {
                currentRow.emplace_back(x, y, 0.0f);
                lastY = y;
            }
            else {
                rows.push_back(currentRow); // 一行结束
                currentRow.clear();
                currentRow.emplace_back(x, y, 0.0f);
                lastY = y;
            }
        }
    }
    if (!currentRow.empty()) rows.push_back(currentRow);

    // 现在rows是所有行：每隔stepY行，每隔stepX个点采样
    for (int rowIdx = 0; rowIdx < (int)rows.size(); rowIdx += stepY) {
        const auto& row = rows[rowIdx];
        for (int colIdx = 0; colIdx < (int)row.size(); colIdx += stepX) {
            pts.push_back(row[colIdx]);
        }
    }

    std::cout << "[Debug] Loaded sparse terrain points: " << pts.size() << std::endl;
    return pts;
}


void rescaleToBox(std::vector<zVector>& pts, float boxHalfSize = 98.0f)
{
    if (pts.empty()) return;

    zVector minPt = pts[0], maxPt = pts[0];
    for (const auto& p : pts) {
        minPt.x = std::min(minPt.x, p.x);
        minPt.y = std::min(minPt.y, p.y);
        minPt.z = std::min(minPt.z, p.z);
        maxPt.x = std::max(maxPt.x, p.x);
        maxPt.y = std::max(maxPt.y, p.y);
        maxPt.z = std::max(maxPt.z, p.z);
    }

    zVector center = (minPt + maxPt) * 0.5f;
    zVector size = maxPt - minPt;
    float maxDim = std::max({ size.x, size.y, size.z });
    float scale = (maxDim > 1e-6f) ? (boxHalfSize / (maxDim * 0.5f)) : 1.0f;

    for (auto& p : pts) {
        p = (p - center) * scale;
    }
}


static void loadYellowLinesFromCSV(
    const std::string& filename,
    std::vector<std::pair<zSpace::zVector, zSpace::zVector>>& out)
{
    out.clear();
    std::ifstream in(filename);
    if (!in.is_open()) {
        std::cerr << "Cannot open " << filename << "\n";
        return;
    }
    std::string line;
    std::vector<zSpace::zVector> current;
    std::vector<std::vector<zSpace::zVector>> ridges;
    while (std::getline(in, line)) {
        if (line.empty()) {
            if (current.size() >= 2) ridges.push_back(current);
            current.clear();
            continue;
        }
        std::stringstream ss(line);
        float x, y, z; char comma;
        ss >> x >> comma >> y >> comma >> z;
        current.emplace_back(x, y, z);
    }
    if (current.size() >= 2) ridges.push_back(current);

    // 在这里做 2.3 倍放大 + 整体 X 轴负方向平移 5
    for (auto& ridge : ridges)
    {
        for (size_t i = 0; i + 1 < ridge.size(); ++i)
        {
            auto a = ridge[i];
            auto b = ridge[i + 1];

            // 缩放
            a.x *= 2.3f;  a.y *= 2.3f;
            b.x *= 2.3f;  b.y *= 2.3f;
            // 平移 X 轴负方向 5
            a.x -= 5.0f; a.y -= 5.0f;
            b.x -= 5.0f; b.y -= 5.0f;

            out.emplace_back(a, b);
        }
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
            zVector(0, -100, cpt.z),  // bottom
            zVector(0, 100, cpt.z),   // top
            zVector(-100, 0, cpt.z),  // left
            zVector(100, 0, cpt.z)    // right
        };

        // 2. 固定 seed1 为两个定值，直接构造 rays 和 seed1
        rays.clear();
        seed1.clear();

        // 第一个种子
        seed1.emplace_back(65.0f, 8.0f, 0.0f);
        rays.emplace_back(cpt, seed1[0]);
        std::cout << "seed0 = (65, 25, 0)\n";

        // 第二个种子
        seed1.emplace_back(-48.0f, 58.0f, 0.0f);
        rays.emplace_back(cpt, seed1[1]);
        std::cout << "seed1 = (-25, 25, 0)\n";

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
        const float minX = -100.f, maxX = 100.f;
        const float minY = -100.f, maxY = 100.f;
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
        ellipses.push_back({ seed0, 10.f, 10.f, theta0 });
        // 四个固定中心点
        for (int i = 0; i < seed1.size(); ++i) {
            const zVector& pt = seed1[i];
            const auto& ray = rays[i];

            zVector dir = ray.second - ray.first;
            float theta = std::atan2(dir.y, dir.x);  // 方向角（xy平面）

            // 添加椭圆，a=12长轴，b=6短轴，theta 为朝向方向
            ellipses.push_back({ pt, 10.f, 10.f, theta });
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



    //—————————————————————————LINE —————————————————————————————//
    //————————————————————————————————————————————————————————//
    //————————————————————————————————————————————————————————//
    //————————————————————————————————————————————————————————//
    //————————————————————————————————————————————————————————//
    //————————————————————————————————————————————————————————//





    //—————————————————————————LINE —————————————————————————————//
    //————————————————————————————————————————————————————————//
    //————————————————————————————————————————————————————————//
    //————————————————————————————————————————————————————————//
    //————————————————————————————————————————————————————————//
    //————————————————————————————————————————————————————————//



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



/* ---------- PointCloudNet : Load & Draw Point Cloud as Curve ---------- */
class PointCloudNet
{
public:
    std::vector<zVector> points;  // Original loaded points
    std::vector<zVector> curvePoints;  // Rescaled curve points
    std::vector<std::pair<zVector, zVector>> segments;

    // Load points from a file (format: x,y,z per line)
    bool loadFromFile(const std::string& filename)
    {
        points.clear();
        curvePoints.clear();

        std::ifstream inFile(filename);
        if (!inFile.is_open())
        {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return false;
        }

        std::string line;
        while (std::getline(inFile, line))
        {
            std::istringstream iss(line);
            std::string xStr, yStr, zStr;

            // Split line by comma
            if (std::getline(iss, xStr, ',') && std::getline(iss, yStr, ',') && std::getline(iss, zStr, ','))
            {
                float x = std::stof(xStr);
                float y = std::stof(yStr);
                float z = std::stof(zStr);
                points.push_back(zVector(x, y, z));
            }
        }

        std::cout << "Loaded " << points.size() << " points from " << filename << std::endl;

        // Automatically compute bounding box + rescale
        computeBoundingBoxAndRescale();

        return true;
    }

    // Compute bounding box and rescale points to [-50,50] cube
    void computeBoundingBoxAndRescale()
    {
        if (points.empty()) return;

        // Compute bounding box
        zVector minPt = points[0];
        zVector maxPt = points[0];

        for (const auto& pt : points)
        {
            minPt.x = std::min(minPt.x, pt.x);
            minPt.y = std::min(minPt.y, pt.y);
            minPt.z = std::min(minPt.z, pt.z);

            maxPt.x = std::max(maxPt.x, pt.x);
            maxPt.y = std::max(maxPt.y, pt.y);
            maxPt.z = std::max(maxPt.z, pt.z);
        }

        // Print bounding box
        std::cout << "Bounding box: min(" << minPt.x << ", " << minPt.y << ", " << minPt.z << ") "
            << "max(" << maxPt.x << ", " << maxPt.y << ", " << maxPt.z << ")" << std::endl;

        // Compute center and scale
        zVector center = (minPt + maxPt) * 0.5f;
        zVector size = maxPt - minPt;

        // Determine max dimension to fit uniformly in [-50,50]
        float maxDim = std::max({ size.x, size.y, size.z });
        float targetHalfSize = 100.0f;
        float scale = (maxDim > 1e-6f) ? (targetHalfSize / (maxDim * 0.5f)) : 1.0f;

        std::cout << "Scaling factor: " << scale << std::endl;

        // Rescale + center points
        curvePoints.clear();
        for (const auto& pt : points)
        {
            zVector shifted = pt - center;
            zVector scaled = shifted * scale;
            curvePoints.push_back(scaled);
        }


    }

    // Draw the rescaled curve
    void drawCurve() const
    {
        if (curvePoints.empty()) return;

        glColor3f(0.2f, 0.7f, 0.2f);  // Green curve
        glLineWidth(3.0f);
        glBegin(GL_LINE_STRIP);
        for (const auto& pt : curvePoints)
        {
            glVertex3f(pt.x, pt.y, pt.z);
        }
        glEnd();
        glLineWidth(1.0f);
    }

    // Optional: Draw points as dots (rescaled points)
    //void drawPoints() const
    //{
    //    if (curvePoints.empty()) return;

    //    glColor3f(1.0f, 0.0f, 0.0f);  // Red points
    //    glPointSize(5.0f);
    //    glBegin(GL_POINTS);
    //    for (const auto& pt : curvePoints)
    //    {
    //        glVertex3f(pt.x, pt.y, pt.z);
    //    }
    //    glEnd();
    //    glPointSize(1.0f);
    //}

    void scaleFromTopLeftCorner(float scaleFactor = 0.869f)
    {
        zVector corner(-100.0f, 100.0f, 0.0f);

        for (auto& pt : curvePoints)
        {
            zVector relative = pt - corner;
            relative *= scaleFactor;
            pt = corner + relative;
        }


        // 2. 重新生成 “segments”
        segments.clear();
        if (curvePoints.size() >= 2)
        {
            segments.reserve(curvePoints.size() - 1);
            for (size_t i = 1; i < curvePoints.size(); ++i)
            {
                segments.emplace_back(curvePoints[i - 1], curvePoints[i]);
            }
        }


        std::cout << "Scaled curvePoints from top-left corner by factor " << scaleFactor << std::endl;
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

   //#include "Road.h"


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
    PointCloudNet pointCloud;
    std::vector<zVector> terrainPoints;
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
        // 1) 加载并缩放点云
        pointCloud.loadFromFile("data/Land Border.txt");
        pointCloud.scaleFromTopLeftCorner(0.896f);
        terrainPoints = loadTerrainPoints("data/02_Terrain Points.txt");
        rescaleToBox(terrainPoints, 100.0f);
        // 2) 闭合折线，让它首尾相连成多边形
        {
            auto& segs = pointCloud.segments;
            if (!segs.empty()) {
                // 最后一条线段的终点连回第一条的起点
                segs.emplace_back(segs.back().second, segs.front().first);
            }
        }

        // 3) 传给 EllipseExpansion 和 ColorGrowth
        ellipses.setPointSegments(pointCloud.segments);
        growth.setBoundarySegments(pointCloud.segments);

        // 4) 其它初始化——随机种子、模式开关
        globalRng.seed(currentSeed);
        showEllipseMode = false;

        loadYellowLinesFromCSV("data/smoothed_main_ridges.csv", yellowLines);
        std::cerr << "[Debug] scene.yellowLines=" << yellowLines.size() << std::endl;

        // 5) 重建网格、交点
        net.build();
        marks.build(net);

        // 6) 重置 growth，并添加主种子
        growth.clear();
        growth.setPlacementMode(ColorGrowth::PlacementMode::Contour);
        growth.addSeed(0, cpt);

        for (auto& z : net.seed1) {
            growth.addSeed(1, z);
        }

        // 8) 生成等高线并注入 growth / ellipses
        net.contour0Lines.clear();
        net.contour1Lines.clear();
        net.drawContourMap(cpt, marks.intPts, 120, 40, 0.45);
        growth.setContour0Lines(net.contour0Lines);
        ellipses.setContour1Lines(net.contour1Lines);

        // 9) 回到初始状态
        stage = 0;
    }




    void next() {
        if (showEllipseMode) return;
        stage = (stage + 1) % 3;
        if (stage == 1) {
            // 只在进入 Stage 1 的时候，把硬写的 seed1 加进来
            for (auto& z : net.seed1) {
                growth.addSeed(1, z);

            }

        }
    }


    void draw() {

        pointCloud.drawCurve();
        // pointCloud.drawPoints();

        glColor3f(0.15f, 0.22f, 0.6f);
        glPointSize(4.0f);
        glBegin(GL_POINTS);
        for (auto& p : terrainPoints) {
            if (!growth.pointInBoundary(p))   // 只画外部
                glVertex3f(p.x, p.y, p.z);
        }
        glEnd();





        // —— 一上来就画平滑主脊线段 —— 
        glColor3f(1.0f, 1.0f, 1.0f);    // 黄色
        glLineWidth(2.0f);
        for (auto& seg : yellowLines) {
            drawLine(z2A(seg.first), z2A(seg.second));
        }
        glLineWidth(1.0f);

        // —— T 模式优先：只画三角形 —— 
        if (showTriangleOnly) {
            ellipses.drawVectorTriangles();
            return;
        }

        if (showEllipseMode) {
            // 1) 先画等高线（如需要）
            // net.drawContourMap(cpt, marks.intPts, 120, 40, 0.35f);
            // 2) 更新并画椭圆
            ellipses.update();
            ellipses.draw();

            // 画每个椭圆的分区中心点
            for (const auto& ep : ellipses.ellipses) {
                for (const auto& c : ep.regionCenters) {
                    glColor3f(1, 0, 0);
                    drawCircle(z2A(c), 1.5f, 16);
                }
            }

            // 其它如有黄线逻辑按需
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

        // 普通模式下的其他绘制
        net.drawContourMap(cpt, marks.intPts, 120, 40, 0.45);
        // net.drawCurves();
        // net.drawRays();
        if (stage == 1) {
            // Stage 1：绘制我们手写的两个 seed1
            net.drawSeeds();

        }
        else {
            // 其它阶段不绘制 marks

        }

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
                // ellipses.splitAllEllipsesWithContourAndVector();
                cmdLog += 'a';
            }
            break;
        case 'w': case 'W':
            if (showEllipseMode) {
                // 只对小椭圆（circles）做一次三点平滑
                ellipses.circleCornerSmooth(SMOOTH_PASSES);
                cmdLog += 'w';
                glutPostRedisplay();
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
            net.drawContourMap(cpt, marks.intPts, 120, 40, 0.45f);
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
            growth.removeNearNetwork(yellowLines, 2.6f);
            growth.removeNearNetwork(pointCloud.segments, 2.6f);

            recordSnapshot();
            break;
        }
        case 'm': case 'M':
            // 删除所有距紫色点半径 12 单位内的异色点
            growth.removePointsNearPurple(6.0f);
            cmdLog += 'm';
            break;

        case 'o': case 'O':
            if (showEllipseMode) {
                ellipses.toggleOffset();
                cmdLog += 'o';
            }
            break;
        case 'q': case 'Q': {
            if (showEllipseMode) {
                // 1) 把每个大椭圆按当前 boundary 切分，计算出 regionCenters
                ellipses.splitAllEllipsesWithContourAndVector();

                // 2) 根据 regionCenters 生成小圆检测种子，并开启扩散（内部会设置 deforming=true）
                ellipses.initCircleSeeds();

                // 3) （可选）记录日志、请求重绘
                cmdLog += 'q';
                glutPostRedisplay();
            }
            break;
        }

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
