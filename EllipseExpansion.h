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
        std::vector<zVector> regionCenters;
        float                    maxRadius;  // 包围半径
    };

    // ———— 椭圆扩散类 ————
    class EllipseExpansion {
    public:
        //—— 私有成员 ——
        std::vector<EPoint>                           ellipses;
        std::vector<std::pair<zVector, zVector>>       splineLines, rayLines, edgeLines;
        std::vector<std::pair<zVector, zVector>>       yellowLines;                           // 新增黄线存储
        std::vector<std::vector<zVector>>             contour0Lines;
        std::vector<std::vector<zVector>>             contour1Lines;
        std::vector<std::pair<zVector, zVector>> pointSegments;

        const int                                     STEPS = 45;
        const float                                   SAME_DIST = 1.2;
        const float                                   DIFF_DIST = 1.2f;
        const std::array<float, 6>                    growRates = { 0.3f,0.3f,0.3f,0.3f,0.3f,0.3f };
        bool                                          showOffset = false;
        bool                                          deforming = false;
        bool                                          fillMode = false;
        double                                        lastTime = 0.0;
        double lastTimeCircles = 0.0;   // 驱动小圆扩散
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


        struct EpDrawData {
            std::vector<std::pair<zVector, zVector>> lines;
            std::vector<zVector> midpoints;
        };
        std::vector<EpDrawData> drawDataPerEllipse;

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
                for (auto& contour : contour0Lines) {
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
                zVector major = tangent * 1.6f;
                zVector minor = zVector(-tangent.y, tangent.x, 0) * 0.8f;

                // 3) 初始化 EPoint
                EPoint ep;
                ep.pos = p;
                ep.rep = minor;
                ep.col = gp.col;
                ep.boundary.resize(STEPS);
                ep.frozen.assign(STEPS, false);
                ep.maxRadius = 20.0f;

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
            contour0Lines = evens;
        }

        void setPointSegments(const std::vector<std::pair<zVector, zVector>>& segs) {
            pointSegments = segs;
        }


        //-----------------------------------------



        bool segmentsIntersect(const zVector& a, const zVector& b, const zVector& c, const zVector& d, zVector& intersection) {
            // 2D平面直线段相交
            float s1_x = b.x - a.x;
            float s1_y = b.y - a.y;
            float s2_x = d.x - c.x;
            float s2_y = d.y - c.y;

            float denom = (-s2_x * s1_y + s1_x * s2_y);
            if (fabs(denom) < 1e-8f) return false; // 平行

            float s = (-s1_y * (a.x - c.x) + s1_x * (a.y - c.y)) / denom;
            float t = (s2_x * (a.y - c.y) - s2_y * (a.x - c.x)) / denom;
            if (s >= 0 && s <= 1 && t >= 0 && t <= 1) {
                intersection.x = a.x + (t * s1_x);
                intersection.y = a.y + (t * s1_y);
                intersection.z = a.z; // 忽略z
                return true;
            }
            return false;
        }


        zVector polygonCentroid(const std::vector<zVector>& poly) {
            double A = 0, Cx = 0, Cy = 0;
            int n = poly.size();
            for (int i = 0; i < n; ++i) {
                const auto& p0 = poly[i];
                const auto& p1 = poly[(i + 1) % n];
                double cross = p0.x * p1.y - p1.x * p0.y;
                A += cross;
                Cx += (p0.x + p1.x) * cross;
                Cy += (p0.y + p1.y) * cross;
            }
            A *= 0.5;
            if (std::abs(A) < 1e-7) return zVector(0, 0, 0);
            Cx /= (6 * A);
            Cy /= (6 * A);
            return zVector(Cx, Cy, 0);
        }



        // 主入口：分割并计算中心点
        void EllipseExpansion::splitAllEllipsesWithContourAndVector() {
            for (auto& ep : ellipses) {
                ep.regionCenters.clear();

                zVector center = ep.pos;
                zVector dir = normalize(ep.rep);           // 主方向
                zVector orth = zVector(-dir.y, dir.x, 0.0f); // 正交方向

                // 1) 找 orth（过 center）与 boundary 的两个交点 A、B
                std::vector<zVector> orthoInters;
                std::vector<int>     orthoIdxs;
                for (int i = 0; i < ep.boundary.size(); ++i) {
                    zVector a = ep.boundary[i];
                    zVector b = ep.boundary[(i + 1) % ep.boundary.size()];
                    zVector o0 = center - orth * 2000.f;
                    zVector o1 = center + orth * 2000.f;
                    zVector interP;
                    if (segmentsIntersect(a, b, o0, o1, interP)) {
                        orthoInters.push_back(interP);
                        orthoIdxs.push_back(i + 1);
                    }
                }
                auto get2closest = [&](std::vector<zVector>& pts, std::vector<int>& idxs, const zVector& ref) {
                    if (pts.size() > 2) {
                        std::vector<std::pair<float, int>> dist_idx;
                        for (int i = 0; i < pts.size(); ++i) {
                            float d2 = (pts[i].x - ref.x) * (pts[i].x - ref.x)
                                + (pts[i].y - ref.y) * (pts[i].y - ref.y);
                            dist_idx.emplace_back(d2, i);
                        }
                        std::sort(dist_idx.begin(), dist_idx.end());
                        pts = { pts[dist_idx[0].second], pts[dist_idx[1].second] };
                        idxs = { idxs[dist_idx[0].second], idxs[dist_idx[1].second] };
                    }
                    };
                get2closest(orthoInters, orthoIdxs, center);
                if (orthoInters.size() != 2) continue;

                // 2) 计算 A、B 的中点 center2
                zVector center2 = (orthoInters[0] + orthoInters[1]) * 0.5f;

                // 3) 找经过 center2、方向为 dir（即正交于 orth）那条线与 boundary 的两个交点 C、D
                std::vector<zVector> dir2Inters;
                std::vector<int>     dir2Idxs;
                for (int i = 0; i < ep.boundary.size(); ++i) {
                    zVector a = ep.boundary[i];
                    zVector b = ep.boundary[(i + 1) % ep.boundary.size()];
                    zVector d0 = center2 - dir * 2000.f;
                    zVector d1 = center2 + dir * 2000.f;
                    zVector interP;
                    if (segmentsIntersect(a, b, d0, d1, interP)) {
                        dir2Inters.push_back(interP);
                        dir2Idxs.push_back(i + 1);
                    }
                }
                get2closest(dir2Inters, dir2Idxs, center2);
                if (dir2Inters.size() != 2) continue;

                // 4) 把这四个交点插入到 boundary，得到 newBoundary 和对应索引 cutIdxs
                std::vector<std::pair<int, zVector>> allInserts;
                for (int i = 0; i < 2; ++i) allInserts.emplace_back(orthoIdxs[i], orthoInters[i]);
                for (int i = 0; i < 2; ++i) allInserts.emplace_back(dir2Idxs[i], dir2Inters[i]);
                std::sort(allInserts.begin(), allInserts.end(),
                    [](auto& a, auto& b) { return a.first < b.first; });

                std::vector<zVector> newBoundary = ep.boundary;
                std::vector<int>     cutIdxs;
                int insertCount = 0;
                for (auto& pr : allInserts) {
                    newBoundary.insert(newBoundary.begin() + pr.first + insertCount,
                        pr.second);
                    cutIdxs.push_back(pr.first + insertCount);
                    ++insertCount;
                }

                // 5) 用四个插入点把 newBoundary 分成 4 块，每块先加 center2，再加弧段顶点，最后算质心
                int sz = newBoundary.size();
                std::sort(cutIdxs.begin(), cutIdxs.end());
                for (int s = 0; s < 4; ++s) {
                    int i0 = cutIdxs[s];
                    int i1 = cutIdxs[(s + 1) % 4];

                    // 构造扇区多边形
                    std::vector<zVector> region;
                    region.push_back(center2);  // 关键：把 center2 作为第一个顶点

                    // 把 boundary 上从 i0 (含) 到 i1 (含) 的顶点依次加入
                    int k = i0;
                    do {
                        region.push_back(newBoundary[k]);
                        k = (k + 1) % sz;
                    } while (k != (i1 + 1) % sz);

                    if (region.size() >= 3) {
                        // 直接求出质心
                        zVector c = polygonCentroid(region);
                        // 如果需要可以再做一次往中心的插值：
                        // c = center + (c - center) * 0.5f;
                        ep.regionCenters.push_back(c);
                    }
                }
            }
        }



        //-----------------------------------------

        //----------------------------------------



        struct CPoint {
            zVector                  pos;       // 圆心
            std::vector<zVector>     boundary;  // 圆周上的检测点
            std::vector<bool>        frozen;    // 每个检测点是否被冻结
            float                    maxRadius; // 最高扩散半径（当前无用，可留着扩展）
        };
        std::vector<CPoint> circles;            // 存放所有质心圆

        // 1) 根据 regionCenters 生成小圆检测种子，并初始化小圆时钟
        inline void EllipseExpansion::initCircleSeeds() {
            circles.clear();
            for (auto& ep : ellipses) {
                for (auto& c : ep.regionCenters) {
                    // 只在 x∈[-50,-10], y∈[10,90] 的正方形区域内生成小圆
                    if (c.x < -50.0f || c.x > -10.0f || c.y < 10.0f || c.y > 90.0f)
                        continue;

                    CPoint cp;
                    cp.pos = c;
                    cp.maxRadius = 0.3f;
                    cp.boundary.resize(STEPS);
                    cp.frozen.assign(STEPS, false);
                    // 在质心周围等分生成检测点
                    for (int i = 0; i < STEPS; ++i) {
                        float θ = 2 * M_PI * i / STEPS;
                        cp.boundary[i] = {
                            c.x + std::cos(θ) * 0.8f,
                            c.y + std::sin(θ) * 0.8f,
                            c.z
                        };
                    }
                    circles.push_back(cp);
                }
            }
            // 使能变形/扩散逻辑
            deforming = true;
            // 重置小圆专用时钟
            lastTimeCircles = std::clock() / double(CLOCKS_PER_SEC);
        }

        inline void EllipseExpansion::updateCircles() {
            if (!deforming) return;
            double now = std::clock() / double(CLOCKS_PER_SEC);
            float dt = float(now - lastTimeCircles);
            lastTimeCircles = now;

            for (auto& cp : circles) {
                for (int k = 0; k < STEPS; ++k) {
                    if (cp.frozen[k]) continue;
                    zVector P = cp.boundary[k];

                    // —— 冻结检测 ——  
                    bool freeze = false;

                    // 2.1 撞大椭圆边
                    for (auto& ep : ellipses) {
                        for (int j = 0; j + 1 < (int)ep.boundary.size(); ++j) {
                            if (distPointSegment(P, ep.boundary[j], ep.boundary[j + 1]) < 0.8f) {
                                freeze = true; break;
                            }
                        }
                        if (freeze) break;
                    }
                    if (freeze) { cp.frozen[k] = true; continue; }

                    // 2.2 圆与圆碰撞
                    for (auto& other : circles) {
                        if (&other == &cp) continue;
                        for (auto& Q : other.boundary) {
                            if ((P - Q).length() < 0.8f) {
                                freeze = true; break;
                            }
                        }
                        if (freeze) break;
                    }
                    if (freeze) { cp.frozen[k] = true; continue; }

                    // 2.3 点云线段检测
                    for (auto& seg : pointSegments) {
                        if (distPointSegment(P, seg.first, seg.second) < DIFF_DIST) {
                            freeze = true; break;
                        }
                    }
                    if (freeze) { cp.frozen[k] = true; continue; }

                    // 2.4 偶数层等高线（contour0）检测，阈值 0.8
                    for (auto& contour : contour1Lines) {
                        for (size_t j = 0; j + 1 < contour.size(); j += 2) {
                            if (distPointSegment(P, contour[j], contour[j + 1]) < 0.8f) {
                                freeze = true;
                                break;
                            }
                        }
                        if (freeze) break;
                    }
                    if (freeze) { cp.frozen[k] = true; continue; }

                    // —— 扩散 ——  
                    zVector dir = P - cp.pos;
                    zVector cand = cp.pos + dir * (1.0f + growRates[int(ColorType::BLUE)] * dt);
                    // 越界就冻结
                    if (cand.x < -100 || cand.x > 100 || cand.y < -100 || cand.y > 100) {
                        cp.frozen[k] = true;
                    }
                    else {
                        cp.boundary[k] = cand;
                    }
                }
            }
        }


        // 3) 小圆的绘制保持不变
        inline void EllipseExpansion::drawCircles() {
            // 3.1 画圆轮廓
            glColor3f(1, 0, 0);
            for (auto& cp : circles) {
                glBegin(GL_LINE_LOOP);
                for (auto& p : cp.boundary) glVertex3f(p.x, p.y, p.z);
                glEnd();
            }
            // 3.2 画检测点（小黑点）
            glColor3f(0, 0, 0);
            for (auto& cp : circles) {
                for (auto& p : cp.boundary) {
                    drawCircle(z2A(p), 0.05f, 6);
                }
            }
        }


        void circleCornerSmooth(int passes)
        {
            const int N = STEPS;
            for (int it = 0; it < passes; ++it)
            {
                for (auto& cp : circles)
                {
                    auto bak = cp.boundary;
                    for (int k = 0; k < N; ++k)
                    {
                        int km = (k - 1 + N) % N;
                        int kp = (k + 1) % N;
                        cp.boundary[k] = bak[k] * 0.75f
                            + bak[km] * 0.125f
                            + bak[kp] * 0.125f;
                    }
                }
            }
        }



        //-------------------------------------------


        //—— 每帧更新 ——
        void update()
        {
            if (!deforming) return;
            double now = std::clock() / double(CLOCKS_PER_SEC);
            float dt = float(now - lastTime);
            lastTime = now;

            for (auto& ep : ellipses) {
                // 1) 免检区域
                float dx0 = ep.pos.x + 100.0f, dy0 = ep.pos.y + 100.0f;
                bool skipDetect = (dx0 * dx0 + dy0 * dy0 <= 0.0f);

                for (int k = 0; k < STEPS; ++k) {
                    if (ep.frozen[k]) continue;
                    zVector oldP = ep.boundary[k];

                    // 2) yellowLines 冻结检测
                    bool freeze = false;
                    if (!skipDetect) {
                        for (auto& seg : yellowLines) {
                            if (distPointSegment(oldP, seg.first, seg.second) < 1.1f) {
                                freeze = true;
                                break;
                            }
                        }
                    }
                    if (freeze) { ep.frozen[k] = true; continue; }

                    //// 3) 偶数层等高线检测
                    //for (auto& contour : contour1Lines) {
                    //    for (size_t j = 0; j + 1 < contour.size(); j += 2) {
                    //        if (distPointSegment(oldP, contour[j], contour[j + 1]) <= 0.1f) {
                    //            freeze = true; break;
                    //     }
                    //    }
                    //    if (freeze) break;
                    //}
                    //if (freeze) { ep.frozen[k] = true; continue; }


                    // ——— 4) 点云线段检测 ———

                    for (auto& seg : pointSegments) {
                        if (distPointSegment(ep.boundary[k], seg.first, seg.second) < 1.5f) {
                            freeze = true;
                            break;
                        }
                    }
                    if (freeze) {
                        ep.frozen[k] = true;
                        continue;
                    }


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
                    if (cand.x < -100 || cand.x > 100 || cand.y < -100 || cand.y > 100) {
                        ep.frozen[k] = true;
                    }
                    else {
                        ep.boundary[k] = cand;
                    }
                }
            }
            updateCircles();

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
                    case RED:    // #87151F → (135, 21, 31)
                        r = 135.0f / 255.0f;  // ≈0.529f
                        g = 21.0f / 255.0f;  // ≈0.082f
                        b = 31.0f / 255.0f;  // ≈0.122f
                        break;
                    case ORANGE:      // #F27D20 → (242, 125, 32)
                        r = 242.0f / 255.0f;  // ≈0.949f
                        g = 100.0f / 255.0f;  // ≈0.490f
                        b = 12.0f / 255.0f;  // ≈0.125f
                        break;
                    case GREEN:        // #80CBA4 → (128, 203, 164)
                        r = 128.0f / 255.0f;  // ≈0.502f
                        g = 203.0f / 255.0f;  // ≈0.796f
                        b = 164.0f / 255.0f;  // ≈0.643f
                        break;
                    case YELLOW:  // #7D CF DB → (125, 207, 219)
                        r = 125.0f / 255.0f;  // ≈0.490f
                        g = 207.0f / 255.0f;  // ≈0.812f
                        b = 219.0f / 255.0f;  // ≈0.859f
                        break;
                    case BLUE:        // #3683AF → (54, 131, 175)
                        r = 54.0f / 255.0f;  // ≈0.212f
                        g = 131.0f / 255.0f;  // ≈0.514f
                        b = 175.0f / 255.0f;  // ≈0.686f
                        break;
                    case PURPLE:        // #28336B → (40, 51, 107)
                        r = 40.0f / 255.0f;  // ≈0.157f
                        g = 51.0f / 255.0f;  // ≈0.200f
                        b = 107.0f / 255.0f;  // ≈0.420f
                        break;
                    /*case BLUE:   r = 0.220f; g = 0.906f; b = 0.890f; break;
                    case RED:    r = 0.557f; g = 0.000f; b = 0.000f; break;
                    case PURPLE: r = 0.506f; g = 0.000f; b = 0.800f; break;
                    case GREEN:  r = 0.500f; g = 1.000f; b = 0.500f; break;
                    case YELLOW: r = 1.000f; g = 0.839f; b = 0.200f; break;
                    case ORANGE: r = 0.980f; g = 0.325f; b = 0.024f; break;*/
                    }
                    glColor3f(r, g, b);
                    glBegin(GL_TRIANGLE_FAN);
                    glVertex3f(ep.pos.x, ep.pos.y, ep.pos.z);
                    for (auto& p : ep.boundary) glVertex3f(p.x, p.y, -0.5);
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

            for (const auto& ep : ellipses) {
                glColor3f(1, 0, 0);
                glPointSize(8.0f);
                glBegin(GL_POINTS);
                // 直接把 regionCenters 里的点画出来
                for (const auto& c : ep.regionCenters) {
                    glVertex3f(c.x, c.y, 1);
                }
                glEnd();
            }




            drawCircles();





        }
    };

} // namespace zSpace
