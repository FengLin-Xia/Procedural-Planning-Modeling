#pragma once

#include <vector>
#include <array>
#include <utility>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <headers/zApp/include/zObjects.h>
#include "ColorGrowth.h"  // ���� ColorType��GPoint

namespace zSpace {

    // �������� �������� ��������
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

    // �������� ��Բ�߽��ṹ ��������
    struct EPoint {
        zVector                  pos;        // Բ��
        zVector                  rep;        // ���������������᷽��
        ColorType                col;        // ��ɫ
        std::vector<zVector>     boundary;   // �߽��
        std::vector<bool>        frozen;     // �����־
        float                    maxRadius;  // ��Χ�뾶
    };

    // �������� ��Բ��ɢ�� ��������
    class EllipseExpansion {
        //���� ˽�г�Ա ����
        std::vector<EPoint>                           ellipses;
        std::vector<std::pair<zVector, zVector>>       splineLines, rayLines, edgeLines;
        std::vector<std::pair<zVector, zVector>>       yellowLines;                           // �������ߴ洢
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
        bool    showVectorTriangle = false;      // �������� t ����ʾ vector+������
        float   vectorTriHeight = 2.f;       // �����εġ��ߡ�
        float   vectorBaseWidth = 0.8f;       // �����εġ��ױ߿�ȡ�
        static inline zVector normalize(const zVector& v) {
            float len = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
            if (len > 1e-6f) {
                return zVector{ v.x / len, v.y / len, v.z / len };
            }
            return v;  // ������ԭ������
        }

        zVector                      seed0_;        // ������������
        std::vector<zVector>         seeds_;        // ������net.build() ����������� seed1


    public:


        //���� ��ʼ�������ã�drawContourMap ������ ����
        void init(
            const zVector& seed0,
            const std::vector<zVector>& seed1,
            const std::vector<GPoint>& growthPts,
            const std::vector<std::pair<zVector, zVector>>& _splineLines,
            const std::vector<std::pair<zVector, zVector>>& _rayLines,
            const std::vector<std::pair<zVector, zVector>>& _edgeLines,
            const std::vector<std::pair<zVector, zVector>>& yellowLinesParam  // ��������
        )
        {
            ellipses.clear();
            splineLines = _splineLines;
            rayLines = _rayLines;
            edgeLines = _edgeLines;
            yellowLines = yellowLinesParam;  // �洢��������
            deforming = false;
            fillMode = false;



            // �����߶ȳ��ݶ� lambda
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

                // 1) �ҵ� nearest contour1Lines ���� p ������߶Σ������ķ���������
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

                // 2) ��������߷�������Բ������
                zVector major = tangent * 0.4f;
                zVector minor = zVector(-tangent.y, tangent.x, 0) * 0.8f;

                // 3) ��ʼ�� EPoint
                EPoint ep;
                ep.pos = p;
                ep.rep = minor;
                ep.col = gp.col;
                ep.boundary.resize(STEPS);
                ep.frozen.assign(STEPS, false);
                ep.maxRadius = 7.0f;

                // 4) ���ݳ��������ɱ߽�
                for (int k = 0; k < STEPS; ++k) {
                    float �� = 2 * M_PI * k / STEPS;
                    ep.boundary[k] = p
                        + major * std::cos(��)
                        + minor * std::sin(��);
                }
                ellipses.push_back(ep);
            }
        }


        void zSpace::EllipseExpansion::toggleVectorTriangle() {
            showVectorTriangle = !showVectorTriangle;
        }

        // ��ʵ���ƣ�����ÿ�� ellipse��������᷽�򣬻���+����������
        void zSpace::EllipseExpansion::drawVectorTriangles() {
            if (!showVectorTriangle) return;

            // ���� ������ֱ�Ӷ��������ų�� ���� 
            std::vector<zVector> repulsionPoints;
            repulsionPoints.emplace_back(-45.f, 15.f, 0.f);
            repulsionPoints.emplace_back(19.2046f, 4.86243f, 0.f);
            repulsionPoints.emplace_back(-8.00135f, -38.4425f, 0.f);

            const float H = vectorTriHeight;
            const float W = vectorBaseWidth;

            for (auto& ep : ellipses) {
                const auto& C = ep.pos;

                // 1) �������������ķ�ƽ������
                zVector rep{ 0,0,0 };
                for (auto& S : repulsionPoints) {
                    zVector d{ C.x - S.x, C.y - S.y, 0 };
                    float d2 = d.x * d.x + d.y * d.y;
                    if (d2 > 1e-6f) {
                        rep.x += d.x / d2;
                        rep.y += d.y / d2;
                    }
                }

                // 2) ��һ����������
                float len = std::sqrt(rep.x * rep.x + rep.y * rep.y);
                if (len < 1e-6f) continue;
                //zVector dir{ -rep.y / len, rep.x / len, 0.0f };
                zVector dir{ rep.x / len, rep.y / len, 0.0f };


                // 3) ���������ζ���͵ױ�
                zVector apex{ C.x + dir.x * H, C.y + dir.y * H, C.z };
                zVector perp{ -dir.y, dir.x, 0 };
                zVector b1{ C.x + perp.x * (W * 0.5f), C.y + perp.y * (W * 0.5f), C.z };
                zVector b2{ C.x - perp.x * (W * 0.5f), C.y - perp.y * (W * 0.5f), C.z };

                // 4) ����������
                glColor3f(0, 0, 0);
                glLineWidth(2.0f);
                glBegin(GL_LINES);
                glVertex3f(C.x, C.y, C.z);
                glVertex3f(apex.x, apex.y, apex.z);
                glEnd();

                // 5) ������������
                glColor3f(0, 0, 0);
                glBegin(GL_TRIANGLES);
                glVertex3f(apex.x, apex.y, apex.z);
                glVertex3f(b1.x, b1.y, b1.z);
                glVertex3f(b2.x, b2.y, b2.z);
                glEnd();
            }
        }




        //���� �л�������ģʽ ����
        void startDeform()
        {
            deforming = true;
            lastTime = std::clock() / double(CLOCKS_PER_SEC);
        }

        //���� ����ƽ���������� ����
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

        //���� ����ż����ȸ������ڶ����� ����
        void setContour1Lines(const std::vector<std::vector<zVector>>& evens)
        {
            contour1Lines = evens;
        }

        //���� ÿ֡���� ����
        void update()
        {
            if (!deforming) return;
            double now = std::clock() / double(CLOCKS_PER_SEC);
            float dt = float(now - lastTime);
            lastTime = now;

            for (auto& ep : ellipses) {
                // 1) �������
                float dx0 = ep.pos.x + 50.0f, dy0 = ep.pos.y + 50.0f;
                bool skipDetect = (dx0 * dx0 + dy0 * dy0 <= 0.0f);

                for (int k = 0; k < STEPS; ++k) {
                    if (ep.frozen[k]) continue;
                    zVector oldP = ep.boundary[k];

                    // 2) yellowLines ������
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

                    // 3) ż����ȸ��߼��
                   // for (auto& contour : contour1Lines) {
                       // for (size_t j = 0; j + 1 < contour.size(); j += 2) {
                           // if (distPointSegment(oldP, contour[j], contour[j + 1]) <= 0.1f) {
                              //  freeze = true; break;
                         // }
                       // }
                        //if (freeze) break;
                  //  }
                   // if (freeze) { ep.frozen[k] = true; continue; }

                    // 4) ��Բ����ײ���
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

                    // 5) ��ɢ�ƶ�
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

        //���� �л���ƫ����ʾ ����
        void toggleOffset()
        {
            showOffset = !showOffset;
        }

        //���� ����������Բ ����
        void draw()
        {
            for (auto& ep : ellipses) {
                if (fillMode) {
                    // ��ɫ���
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

                    // ��ƫ�����
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

                // �������������
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
                // ֻ�������Σ������������л���
                drawVectorTriangles();
                return;
            }
        }
    };

} // namespace zSpace