#pragma once

#include <vector>
#include <array>
#include <utility>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <headers/zApp/include/zObjects.h>
#include "ColorGrowth.h"  // 定义 ColorType、GPoint

namespace zSpace {

    // ———— 辅助函数 ————
    inline float distPointSegment(
        const zVector& P,
        const zVector& A,
        const zVector& B)
    {
        zVector AB = B - A;
        float len2 = AB.x * AB.x + AB.y * AB.y;
        if (len2 < 1e-6f) return (P - A).length();
        float t = ((P.x - A.x) * AB.x + (P.y - A.y) * AB.y) / len2;
        t = std::fmax(0.0f, std::fmin(1.0f, t));
        zVector H = A + AB * t;
        return (P - H).length();
    }

    // ———— 椭圆边界点结构 ————
    struct EPoint {
        zVector                  pos;        // 圆心
        zVector                  rep;        // 主方向向量（短轴方向）
        ColorType                col;        // 颜色
        std::vector<zVector>     boundary;   // 边界点
        std::vector<bool>        frozen;     // 冻结标志
        float                    maxRadius;  // 包围半径
    };

    // ———— 椭圆扩散类 ————
    class EllipseExpansion {
        //—— 私有成员 ——
        std::vector<EPoint>                           ellipses;
        std::vector<std::pair<zVector, zVector>>       splineLines, rayLines, edgeLines;
        std::vector<std::pair<zVector, zVector>>       yellowLines;                           // 新增黄线存储
        std::vector<std::vector<zVector>>             contour1Lines;
        const int                                     STEPS = 35;
        const float                                   SAME_DIST = 0.8;
        const float                                   DIFF_DIST = 0.8f;
        const std::array<float, 6>                    growRates = { 0.3f,0.3f,0.3f,0.3f,0.3f,0.3f };
        bool                                          showOffset = false;
        bool                                          deforming = false;
        bool                                          fillMode = false;
        double                                        lastTime = 0.0;
        std::vector<std::vector<zVector>>             evenContourLines;
        bool    showVectorTriangle = false;      // 新增：按 t 键显示 vector+三角形
        float   vectorTriHeight = 2.f;       // 三角形的“高”
        float   vectorBaseWidth = 0.8f;       // 三角形的“底边宽度”
        static inline zVector normalize(const zVector& v) {
            float len = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
            if (len > 1e-6f) {
                return zVector{ v.x / len, v.y / len, v.z / len };
            }
            return v;  // 零向量原样返回
        }

        zVector                      seed0_;        // 新增：主种子
        std::vector<zVector>         seeds_;        // 新增：net.build() 里算的那两个 seed1


    public:


        //—— 初始化并重置（drawContourMap 触发） ——
        void init(
            const zVector& seed0,
            const std::vector<zVector>& seed1,
            const std::vector<GPoint>& growthPts,
            const std::vector<std::pair<zVector, zVector>>& _splineLines,
            const std::vector<std::pair<zVector, zVector>>& _rayLines,
            const std::vector<std::pair<zVector, zVector>>& _edgeLines,
            const std::vector<std::pair<zVector, zVector>>& yellowLinesParam  // 新增参数
        )
        {
            ellipses.clear();
            splineLines = _splineLines;
            rayLines = _rayLines;
            edgeLines = _edgeLines;
            yellowLines = yellowLinesParam;  // 存储黄线数据
            deforming = false;
            fillMode = false;



            // 构建高度场梯度 lambda
            constexpr float sigma2 = 15.0f * 15.0f;
            auto gradH = [&](float x, float y) {
                zVector g{ 0,0,0 };
                auto add = [&](const zVector& s) {
                    float dx = x - s.x, dy = y - s.y;
                    float w = std::exp(-(dx * dx + dy * dy) / (2 * sigma2));
                    g.x += -dx * w;
                    g.y += -dy * w;
                    };
                add(seed0);
                for (auto& s : seed1) add(s);
                return g;
                };

            for (auto& gp : growthPts) {
                zVector p(gp.pos.x, gp.pos.y, gp.pos.z);

                // 1) 找到 nearest contour1Lines 里与 p 最近的线段，以它的方向做切线
                zVector tangent{ 1,0,0 };
                float bestD = 1e6f;
                for (auto& contour : contour1Lines) {
                    for (size_t i = 0; i + 1 < contour.size(); i += 2) {
                        float d = distPointSegment(p, contour[i], contour[i + 1]);
                        if (d < bestD) {
                            bestD = d;
                            tangent = contour[i + 1] - contour[i];
                        }
                    }
                }
                float tl = std::hypot(tangent.x, tangent.y);
                if (tl > 1e-6f) tangent.x /= tl, tangent.y /= tl;

                // 2) 用这个切线方向构造椭圆主副轴
                zVector major = tangent * 0.4f;
                zVector minor = zVector(-tangent.y, tangent.x, 0) * 0.8f;

                // 3) 初始化 EPoint
                EPoint ep;
                ep.pos = p;
                ep.rep = minor;
                ep.col = gp.col;
                ep.boundary.resize(STEPS);
                ep.frozen.assign(STEPS, false);
                ep.maxRadius = 7.0f;

                // 4) 根据长短轴生成边界
                for (int k = 0; k < STEPS; ++k) {
                    float θ = 2 * M_PI * k / STEPS;
                    ep.boundary[k] = p
                        + major * std::cos(θ)
                        + minor * std::sin(θ);
                }
                ellipses.push_back(ep);
            }
        }


        void zSpace::EllipseExpansion::toggleVectorTriangle() {
            showVectorTriangle = !showVectorTriangle;
        }

        // 其实绘制：遍历每个 ellipse，算出长轴方向，画线+等腰三角形
        void zSpace::EllipseExpansion::drawVectorTriangles() {
            if (!showVectorTriangle) return;

            // —— 在这里直接定义三个排斥点 —— 
            std::vector<zVector> repulsionPoints;
            repulsionPoints.emplace_back(-45.f, 15.f, 0.f);
            repulsionPoints.emplace_back(19.2046f, 4.86243f, 0.f);
            repulsionPoints.emplace_back(-8.00135f, -38.4425f, 0.f);

            const float H = vectorTriHeight;
            const float W = vectorBaseWidth;

            for (auto& ep : ellipses) {
                const auto& C = ep.pos;

                // 1) 计算对这三个点的反平方合力
                zVector rep{ 0,0,0 };
                for (auto& S : repulsionPoints) {
                    zVector d{ C.x - S.x, C.y - S.y, 0 };
                    float d2 = d.x * d.x + d.y * d.y;
                    if (d2 > 1e-6f) {
                        rep.x += d.x / d2;
                        rep.y += d.y / d2;
                    }
                }

                // 2) 归一化合力方向
                float len = std::sqrt(rep.x * rep.x + rep.y * rep.y);
                if (len < 1e-6f) continue;
                //zVector dir{ -rep.y / len, rep.x / len, 0.0f };
                zVector dir{ rep.x / len, rep.y / len, 0.0f };


                // 3) 计算三角形顶点和底边
                zVector apex{ C.x + dir.x * H, C.y + dir.y * H, C.z };
                zVector perp{ -dir.y, dir.x, 0 };
                zVector b1{ C.x + perp.x * (W * 0.5f), C.y + perp.y * (W * 0.5f), C.z };
                zVector b2{ C.x - perp.x * (W * 0.5f), C.y - perp.y * (W * 0.5f), C.z };

                // 4) 画合力向量
                glColor3f(0, 0, 0);
                glLineWidth(2.0f);
                glBegin(GL_LINES);
                glVertex3f(C.x, C.y, C.z);
                glVertex3f(apex.x, apex.y, apex.z);
                glEnd();

                // 5) 画等腰三角形
                glColor3f(0, 0, 0);
                glBegin(GL_TRIANGLES);
                glVertex3f(apex.x, apex.y, apex.z);
                glVertex3f(b1.x, b1.y, b1.z);
                glVertex3f(b2.x, b2.y, b2.z);
                glEnd();
            }
        }




        //—— 切换到变形模式 ——
        void startDeform()
        {
            deforming = true;
            lastTime = std::clock() / double(CLOCKS_PER_SEC);
        }

        //—— 三点平滑，无收缩 ——
        void cornerSmooth(int passes)
        {
            const int N = STEPS;
            for (int it = 0; it < passes; ++it) {
                for (auto& ep : ellipses) {
                    auto bak = ep.boundary;
                    for (int k = 0; k < N; ++k) {
                        int km = (k - 1 + N) % N, kp = (k + 1) % N;
                        ep.boundary[k] = bak[k] * 0.75f
                            + bak[km] * 0.125f
                            + bak[kp] * 0.125f;
                    }
                }
            }
            fillMode = true;
        }

        //—— 传入偶数层等高线用于冻结检测 ——
        void setContour1Lines(const std::vector<std::vector<zVector>>& evens)
        {
            contour1Lines = evens;
        }

        //—— 每帧更新 ——
        void update()
        {
            if (!deforming) return;
            double now = std::clock() / double(CLOCKS_PER_SEC);
            float dt = float(now - lastTime);
            lastTime = now;

            for (auto& ep : ellipses) {
                // 1) 免检区域
                float dx0 = ep.pos.x + 50.0f, dy0 = ep.pos.y + 50.0f;
                bool skipDetect = (dx0 * dx0 + dy0 * dy0 <= 0.0f);

                for (int k = 0; k < STEPS; ++k) {
                    if (ep.frozen[k]) continue;
                    zVector oldP = ep.boundary[k];

                    // 2) yellowLines 冻结检测
                    bool freeze = false;
                    if (!skipDetect) {
                        for (auto& seg : yellowLines) {
                            if (distPointSegment(oldP, seg.first, seg.second) < 1.5f) {
                                freeze = true;
                                break;
                            }
                        }
                    }
                    if (freeze) { ep.frozen[k] = true; continue; }

                    // 3) 偶数层等高线检测
                   // for (auto& contour : contour1Lines) {
                       // for (size_t j = 0; j + 1 < contour.size(); j += 2) {
                           // if (distPointSegment(oldP, contour[j], contour[j + 1]) <= 0.1f) {
                              //  freeze = true; break;
                         // }
                       // }
                        //if (freeze) break;
                  //  }
                   // if (freeze) { ep.frozen[k] = true; continue; }

                    // 4) 椭圆间碰撞检测
                    bool stop = false;
                    for (auto& o : ellipses) {
                        if (&o == &ep) continue;
                        float thresh = (o.col == ep.col ? SAME_DIST : DIFF_DIST);
                        float cd = (ep.pos - o.pos).length();
                        if (cd > ep.maxRadius + o.maxRadius + thresh) continue;
                        for (auto& r : o.boundary) {
                            if ((oldP - r).length() < thresh) { stop = true; break; }
                        }
                        if (stop) break;
                    }
                    if (stop) { ep.frozen[k] = true; continue; }

                    // 5) 扩散移动
                    zVector dir = oldP - ep.pos;
                    zVector cand = ep.pos + dir * (1.0f + growRates[int(ep.col)] * dt);
                    if (cand.x < -50 || cand.x > 50 || cand.y < -50 || cand.y > 50) {
                        ep.frozen[k] = true;
                    }
                    else {
                        ep.boundary[k] = cand;
                    }
                }
            }
        }

        //—— 切换内偏移显示 ——
        void toggleOffset()
        {
            showOffset = !showOffset;
        }

        //—— 绘制所有椭圆 ——
        void draw()
        {
            for (auto& ep : ellipses) {
                if (fillMode) {
                    // 彩色填充
                    float r, g, b;
                    switch (ep.col) {
                    case BLUE:   r = 0.220f; g = 0.906f; b = 0.890f; break;
                    case RED:    r = 0.557f; g = 0.000f; b = 0.000f; break;
                    case PURPLE: r = 0.506f; g = 0.000f; b = 0.800f; break;
                    case GREEN:  r = 0.000f; g = 1.000f; b = 0.000f; break;
                    case YELLOW: r = 1.000f; g = 0.839f; b = 0.200f; break;
                    case ORANGE: r = 0.980f; g = 0.325f; b = 0.024f; break;
                    }
                    glColor3f(r, g, b);
                    glBegin(GL_TRIANGLE_FAN);
                    glVertex3f(ep.pos.x, ep.pos.y, ep.pos.z);
                    for (auto& p : ep.boundary) glVertex3f(p.x, p.y, p.z);
                    glVertex3f(ep.boundary[0].x, ep.boundary[0].y, ep.boundary[0].z);
                    glEnd();

                    // 内偏移填黑
                    if (showOffset) {
                        std::vector<zVector> innerPts;
                        innerPts.reserve(ep.boundary.size());
                        for (auto& p : ep.boundary) {
                            zVector dir = p - ep.pos;
                            float len = std::sqrt(dir.x * dir.x + dir.y * dir.y);
                            if (len > 1e-6f) dir *= (1.0f / len);
                            innerPts.push_back(p - dir * 0.2f);
                        }
                        glColor3f(0, 0, 0);
                        glBegin(GL_TRIANGLE_FAN);
                        glVertex3f(ep.pos.x, ep.pos.y, ep.pos.z);
                        for (auto& ip : innerPts) glVertex3f(ip.x, ip.y, ip.z);
                        glVertex3f(innerPts[0].x, innerPts[0].y, innerPts[0].z);
                        glEnd();
                    }
                }

                // 仅画轮廓或检测点
                if (!showOffset) {
                    glColor3f(1, 1, 1);
                    glBegin(GL_LINE_LOOP);
                    for (auto& p : ep.boundary) glVertex3f(p.x, p.y, p.z);
                    glEnd();
                }
                glColor3f(0, 0, 0);
                for (auto& p : ep.boundary) drawCircle(z2A(p), 0.1f, 8);
            }

            if (showVectorTriangle) {
                // 只画三角形，跳过其它所有绘制
                drawVectorTriangles();
                return;
            }
        }
    };

} // namespace zSpace