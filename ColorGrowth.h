#pragma once

#include <vector>
#include <array>
#include <random>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <limits>
#include <headers/zApp/include/zObjects.h>  // for zVector
#include <headers/zApp/include/zViewer.h>  // for Alice::vec

using namespace zSpace;

// —— 颜色枚举 ——
enum ColorType { BLUE, RED, PURPLE, GREEN, YELLOW, ORANGE };

// —— 生长点结构 ——
struct GPoint {
    Alice::vec pos;
    ColorType col;
    bool active;
    GPoint(Alice::vec p, ColorType c) : pos(p), col(c), active(true) {}
};

// —— 生长连线结构 ——
struct GLine {
    Alice::vec a, b;
    GLine(Alice::vec s, Alice::vec e) : a(s), b(e) {}
};

// —— ColorGrowth 类，支持四种放置模式 ——
class ColorGrowth {
    static constexpr int N_GROUPS = 3;

public:
    // 放置模式
    enum class PlacementMode { ParentRadial, SquareGrid, HexGrid, Contour };

    // 切换模式
    void setPlacementMode(PlacementMode m) { mode = m; }

    // 设置等高线，供 Contour 模式使用
    void setContour0Lines(const std::vector<std::vector<zVector>>& lines_) {
        contour0.clear();
        for (auto& ln : lines_) {
            std::vector<Alice::vec> a;
            a.reserve(ln.size());
            for (auto& p : ln) a.emplace_back(p.x, p.y, p.z);
            contour0.push_back(std::move(a));
        }
    }

    // 设置圆交点，供圆交模式使用
    void setCircleIntersectionPoints(const std::vector<Alice::vec>& pts_) {
        circleIntersectionPoints = pts_;
        usedCirclePoint.assign(pts_.size(), false);
        nextCircleIndex = 0;
    }

    // 恢复／访问状态
    const std::vector<GLine>& getLines() const { return lines; }
    const std::vector<GPoint>& getPoints() const { return pts; }
    void loadState(const std::vector<GPoint>& pts_, const std::vector<GLine>& lines_) {
        pts = pts_;
        lines = lines_;
    }

    // 管理种子和计数
    void clear() {
        for (auto& g : seedsGroups) g.clear();
        pts.clear();
        lines.clear();
        maxCount = 0;
    }
    void addSeed(const zVector& p) { addSeed(0, p); }
    void addSeed(int group, const zVector& p) {
        Alice::vec av(p.x, p.y, p.z);
        seedsGroups[group].push_back(av);
        pts.emplace_back(av, chooseColor(av, BLUE));
    }
    void expandTo50() { maxCount = 150; }
    void expandFull() { maxCount = SIZE_MAX; }

    // 辅助调整
    void recolorNear(const std::vector<Alice::vec>& centers, float radius) {
        float r2 = radius * radius;
        for (auto& c : centers)
            for (auto& gp : pts)
                if ((gp.pos.x - c.x) * (gp.pos.x - c.x) + (gp.pos.y - c.y) * (gp.pos.y - c.y) <= r2)
                    gp.col = BLUE;
    }
    void removeNearNetwork(const std::vector<GLine>& lines_, float radius) {
        float r2 = radius * radius;
        Alice::vec initial = seedsGroups[0].empty() ? Alice::vec(0, 0, 0) : seedsGroups[0][0];
       
        for (int i = int(pts.size()) - 1; i >= 0; --i) {
            auto& gp = pts[i];
            float dx0 = gp.pos.x - initial.x;
            float dy0 = gp.pos.y - initial.y;
           
            bool erase = false;
            for (auto& ln : lines_) {
                Alice::vec a = ln.a, b = ln.b, p = gp.pos;
                Alice::vec ab = b - a, ap = p - a;
                float len2 = ab.x * ab.x + ab.y * ab.y;
                if (len2 < 1e-6f) continue;
                float t = (ab.x * ap.x + ab.y * ap.y) / len2;
                t = t < 0 ? 0 : (t > 1 ? 1 : t);
                Alice::vec proj{ a.x + ab.x * t,a.y + ab.y * t,0 };
                float ddx = p.x - proj.x, ddy = p.y - proj.y;
                if (ddx * ddx + ddy * ddy <= r2) { erase = true; break; }
            }
            if (erase) pts.erase(pts.begin() + i);
        }
    }


    // 新增：删除所有距任意紫色点 radius 范围内的非紫色点
    void removePointsNearPurple(float radius = 15.0f) {
        float r2 = radius * radius;
        // 收集所有紫色点
        std::vector<Alice::vec> purplePos;
        for (auto& gp : pts) {
            if (gp.col == PURPLE) {
                purplePos.push_back(gp.pos);
            }
        }
        // 倒序遍历，遇到非紫色且与任一紫色点距离平方 ≤ r2 则删除
        for (int i = int(pts.size()) - 1; i >= 0; --i) {
            if (pts[i].col == PURPLE) continue;
            bool tooClose = false;
            for (auto& pp : purplePos) {
                float dx = pp.x - pts[i].pos.x;
                float dy = pp.y - pts[i].pos.y;
                float dz = pp.z - pts[i].pos.z;
                if (dx * dx + dy * dy + dz * dz <= r2) {
                    tooClose = true;
                    break;
                }
            }
            if (tooClose) {
                pts.erase(pts.begin() + i);
            }
        }
    }
    void removePointsNearYellow(float radius = 12.0f) {
        float r2 = radius * radius;
        // 收集所有黄色点
        std::vector<Alice::vec> yellowPos;
        for (auto& gp : pts) {
            if (gp.col == YELLOW) {
                yellowPos.push_back(gp.pos);
            }
        }
        // 倒序遍历，遇到非黄色且与任一黄色点距离平方 ≤ r2 则删除
        for (int i = int(pts.size()) - 1; i >= 0; --i) {
            if (pts[i].col == YELLOW) continue;
            bool tooClose = false;
            for (auto& yp : yellowPos) {
                float dx = yp.x - pts[i].pos.x;
                float dy = yp.y - pts[i].pos.y;
                float dz = yp.z - pts[i].pos.z;
                if (dx * dx + dy * dy + dz * dz <= r2) {
                    tooClose = true;
                    break;
                }
            }
            if (tooClose) {
                pts.erase(pts.begin() + i);
            }
        }
    }

    // 统一更新接口
    void update() {
        switch (mode) {
        case PlacementMode::ParentRadial: placeParentRadial(); break;
        case PlacementMode::SquareGrid:  placeOnSquareGrid(); break;
        case PlacementMode::HexGrid:     placeOnHexGrid(); break;
        case PlacementMode::Contour:     placeOnContour();    break;
        }
    }

    // 绘制函数（保留不变）
    void draw() const {
        glPointSize(8.f);
        glBegin(GL_POINTS);
        for (auto& p : pts) {
            setCol(p.col);
            glVertex3f(p.pos.x, p.pos.y, 0);
        }
        glEnd();
        glColor3f(1., 1., 1.);
        glBegin(GL_LINES);
        for (auto& l : lines) {
            glVertex3f(l.a.x, l.a.y, 0);
            glVertex3f(l.b.x, l.b.y, 0);
        }
        glEnd();
    }

private:
    // 当前模式
    PlacementMode mode = PlacementMode::Contour;
    // 核心数据
    std::array<std::vector<Alice::vec>, 3> seedsGroups;
    std::vector<GPoint> pts;
    std::vector<GLine>  lines;
    size_t maxCount = 10;
    std::vector<std::vector<Alice::vec>> contour0;
    std::vector<Alice::vec> circleIntersectionPoints;
    std::vector<bool> usedCirclePoint;
    size_t nextCircleIndex = 0;
    inline static std::mt19937 globalRng{ std::random_device{}() };

    // 工具函数
    static void setCol(ColorType c) {
        switch (c) {
        case BLUE:   glColor3f(0.22f, 0.906f, 0.89f); break;
        case RED:    glColor3f(0.557f, 0, 0);        break;
        case PURPLE: glColor3f(0.506f, 0, 0.8f);    break;
        case GREEN:  glColor3f(0, 0.5f, 0);         break;
        case YELLOW: glColor3f(1, 0.839f, 0.2f);    break;
        case ORANGE: glColor3f(0.98f, 0.325f, 0.024f); break;
        }
    }
    static float dist(const Alice::vec& a, const Alice::vec& b) {
        float dx = a.x - b.x, dy = a.y - b.y;
        return std::sqrt(dx * dx + dy * dy);
    }
    bool within(const Alice::vec& p) const {
        return p.x >= -50 && p.x <= 50 && p.y >= -50 && p.y <= 50;
    }
    bool nearExisting(const Alice::vec& p) const {
        for (auto& q : pts) if (dist(q.pos, p) < 2.5f) return true;
        return false;
    }

    // 颜色继承逻辑
    ColorType chooseColor(const Alice::vec& pos, ColorType parentCol) const {
        std::uniform_real_distribution<float> inheritDist(0.0f, 1.0f);
        if (inheritDist(globalRng) < 0.3f) return parentCol;

        std::array<float, N_GROUPS> minD;
        minD.fill(FLT_MAX);
        for (int g = 0; g < N_GROUPS; ++g) {
            for (auto& s : seedsGroups[g]) {
                minD[g] = std::fmin(minD[g], dist(pos, s));
            }
        }
        if (minD[0] <= 13.0f || minD[1] <= 7.0f || minD[2] <= 2.0f) return BLUE;

        float d = std::sqrt(pos.x * pos.x + pos.y * pos.y);
        std::uniform_int_distribution<int> D(1, 100);
        int r = D(globalRng);

        if (d <= 20.0f) {
            switch (parentCol) {
            case ORANGE:
                if (r <= 60)         return ORANGE;  // 60% → orange
                else if (r <= 80)    return YELLOW;  // 20% → yellow (was 30%)
                else                 return RED;     // 20%
            case YELLOW:
                if (r <= 60)         return ORANGE;  // 60%
                else if (r <= 80)    return YELLOW;  // 20%
                else                 return RED;     // 20%
            case RED:
                if (r <= 70)         return ORANGE;  // 70%
                else if (r <= 90)    return YELLOW;  // 20%
                else                 return RED;     // 10%
            default:
                if (r <= 60)         return ORANGE;
                else if (r <= 80)    return YELLOW;
                else                 return RED;
            }
        }
        else if (d <= 45.0f) {
            switch (parentCol) {
            case ORANGE:
                if (r <= 50)         return ORANGE;  // 50%
                else if (r <= 70)    return YELLOW;  // 20%
                else if (r <= 80)    return RED;     // 10%
                else                 return GREEN;   // 20%
            case YELLOW:
                if (r <= 60)         return ORANGE;  // 60%
                else if (r <= 70)    return YELLOW;  // 10%
                else if (r <= 90)    return GREEN;   // 20%
                else                 return RED;     // 10%
            case RED:
                if (r <= 60)         return ORANGE;  // 60%
                else if (r <= 80)    return YELLOW;  // 20%
                else                 return GREEN;   // 20%
            case GREEN:
                if (r <= 60)         return ORANGE;  // 60%
                else if (r <= 80)    return YELLOW;  // 20%
                else                 return GREEN;   // 20%
            default:
                if (r <= 50)         return ORANGE;
                else if (r <= 70)    return YELLOW;
                else if (r <= 80)    return RED;
                else                 return GREEN;
            }
        }
        else {
            switch (parentCol) {
            case ORANGE:
                if (r <= 70)         return ORANGE;  // 70%
                else if (r <= 75)    return RED;     // 5%
                else if (r <= 95)    return GREEN;   // 20%
                else                 return PURPLE;  // 5%
            case YELLOW:
                if (r <= 70)         return ORANGE;  // 70%
                else if (r <= 80)    return PURPLE;  // 10%
                else if (r <= 90)    return RED;     // 10%
                else                 return GREEN;   // 10%
            case RED:
                if (r <= 70)         return ORANGE;  // 70%
                else if (r <= 80)    return YELLOW;  // 10%
                else                 return GREEN;   // 20%
            case GREEN:
                if (r <= 60)         return ORANGE;  // 60%
                else if (r <= 80)    return YELLOW;  // 20%
                else                 return GREEN;   // 20%
            case PURPLE:
                if (r <= 80)         return PURPLE;  // 80%
                else                 return ORANGE;  // 20%
            default:
                if (r <= 50)         return ORANGE;
                else if (r <= 70)    return YELLOW;
                else if (r <= 75)    return RED;
                else if (r <= 95)    return GREEN;
                else                 return PURPLE;
            }
        }
    }
    // 模式1：父点半径8内随机放置，子点间距>=3
    
    void placeParentRadial() {
        if (pts.size() >= maxCount) return;

        const float parentR = 8.0f, minD = 4.0f;
        std::uniform_real_distribution<float> angDist(0, 2 * M_PI);
        std::uniform_real_distribution<float> unitDist(0, 1);

        std::vector<GPoint> newPts;
        std::vector<GLine>  newLns;

        for (auto& gp : pts) {
            if (!gp.active) continue;

            bool generated = false;
            for (int tries = 0; tries < 10; ++tries) {
                // 均匀圆盘采样
                float ang = angDist(globalRng);
                float r = std::sqrt(unitDist(globalRng)) * parentR;
                Alice::vec c{
                    gp.pos.x + r * std::cos(ang),
                    gp.pos.y + r * std::sin(ang),
                    gp.pos.z
                };
                if (!within(c)) continue;

                // 距离检查：不仅要和旧的 pts 保持距离，也要和本轮 newPts 保持距离
                bool ok = true;
                for (auto& q : pts) {
                    if (dist(q.pos, c) < minD) { ok = false; break; }
                }
                if (!ok) continue;
                for (auto& q : newPts) {
                    if (dist(q.pos, c) < minD) { ok = false; break; }
                }
                if (!ok) continue;

                // 都通过了才算 new point
                newPts.emplace_back(c, chooseColor(c, gp.col));
                newLns.emplace_back(gp.pos, c);
                generated = true;
                break;
            }

            if (!generated) {
                gp.active = false;
            }
        }

        // 推入新点
        for (size_t i = 0; i < newPts.size() && pts.size() < maxCount; ++i) {
            pts.push_back(std::move(newPts[i]));
            lines.push_back(std::move(newLns[i]));
        }
    }



    // 模式2：20×20方格中心放置
    void placeOnSquareGrid() {
        if (pts.size() >= maxCount) return;

        const float half = 50.0f;
        const int   N = 18;                  // 改为 25×25
        const float cell = 100.0f / float(N);    // 每格大小
        const float parentR = 10.0f;

        std::vector<bool> occ(N * N, false);
        std::vector<GPoint> newPts;
        std::vector<GLine>  newLns;

        for (auto& gp : pts) {
            if (!gp.active) continue;

            bool generated = false;
            for (int i = 0; i < N && !generated; ++i) {
                for (int j = 0; j < N && !generated; ++j) {
                    int idx = i * N + j;
                    if (occ[idx]) continue;

                    Alice::vec c{
                        -half + (i + 0.5f) * cell,
                        -half + (j + 0.5f) * cell,
                        gp.pos.z
                    };
                    if (dist(gp.pos, c) > parentR) continue;

                    bool ok = true;
                    for (auto& q : pts) {
                        if (dist(q.pos, c) < cell * 0.5f) { ok = false; break; }
                    }
                    if (!ok) continue;

                    occ[idx] = true;
                    newPts.emplace_back(c, chooseColor(c, gp.col));
                    newLns.emplace_back(gp.pos, c);
                    generated = true;
                }
            }
            if (!generated) gp.active = false;
        }

        for (size_t i = 0; i < newPts.size() && pts.size() < maxCount; ++i) {
            pts.push_back(std::move(newPts[i]));
            lines.push_back(std::move(newLns[i]));
        }
    }



    // 模式3：六边形蜂窝网格中心
    void placeOnHexGrid() {
        if (pts.size() >= maxCount) return;

        const float half = 50, cell = 100.0f / 25.0f;
        const float h = cell * std::sqrt(3.0f) / 2.0f;
        const int rows = int(100.0f / h) + 1;
        const int cols = int(100.0f / cell) + 1;
        const float parentR = 8.0f;

        std::vector<bool> occ(rows * cols, false);
        std::vector<GPoint> newPts;
        std::vector<GLine>  newLns;

        for (auto& gp : pts) {
            if (!gp.active) continue;

            bool generated = false;
            for (int r = 0; r < rows && !generated; ++r) {
                for (int c = 0; c < cols && !generated; ++c) {
                    int idx = r * cols + c;
                    if (occ[idx]) continue;

                    Alice::vec pt{
                        -half + c * cell + ((r & 1) ? cell * 0.5f : 0),
                        -half + r * h,
                        gp.pos.z
                    };
                    if (!within(pt) || dist(gp.pos, pt) > parentR) continue;

                    bool ok = true;
                    for (auto& q : pts) {
                        if (dist(q.pos, pt) < cell * 0.5f) { ok = false; break; }
                    }
                    if (!ok) continue;

                    occ[idx] = true;
                    newPts.emplace_back(pt, chooseColor(pt, gp.col));
                    newLns.emplace_back(gp.pos, pt);
                    generated = true;
                }
            }
            if (!generated) gp.active = false;
        }

        for (size_t i = 0; i < newPts.size() && pts.size() < maxCount; ++i) {
            pts.push_back(std::move(newPts[i]));
            lines.push_back(std::move(newLns[i]));
        }
    }


    // 模式4：原有等高线放置逻辑
    void placeOnContour() {
        if (pts.size() >= maxCount) return;

        std::vector<GPoint> newPts;
        std::vector<GLine>  newLns;
        const float parentRadius = 8.0f;
        const float minDist2 = 16.f;  // 最小距离的平方（3 单位距离）
        size_t orig = pts.size();

        for (size_t i = 0; i < orig && pts.size() + newPts.size() < maxCount; ++i) {
            auto& gp = pts[i];
            if (!gp.active) continue;

            bool generated = false;
            for (auto& line : contour0) {
                for (auto& p : line) {
                    float dx = p.x - gp.pos.x, dy = p.y - gp.pos.y;
                    if (dx * dx + dy * dy > parentRadius * parentRadius) continue;

                    bool ok = true;
                    // 1) 检查与已有点的距离
                    for (auto& q : pts) {
                        float ddx = q.pos.x - p.x, ddy = q.pos.y - p.y;
                        if (ddx * ddx + ddy * ddy < minDist2) {
                            ok = false;
                            break;
                        }
                    }
                    if (!ok) continue;

                    // 2) **新增**：检查与本次循环新生成点的距离
                    for (auto& np : newPts) {
                        float ddx = np.pos.x - p.x, ddy = np.pos.y - p.y;
                        if (ddx * ddx + ddy * ddy < minDist2) {
                            ok = false;
                            break;
                        }
                    }
                    if (!ok) continue;

                    // 满足所有距离约束，生成新点
                    Alice::vec av{ p.x, p.y, p.z };
                    newPts.emplace_back(av, chooseColor(av, gp.col));
                    newLns.emplace_back(gp.pos, av);
                    generated = true;
                    goto DONE;
                }
            }
        DONE:
            if (!generated) gp.active = false;
        }

        // 把 newPts 和 newLns 推入主列表
        for (size_t i = 0; i < newPts.size() && pts.size() < maxCount; ++i) {
            pts.push_back(std::move(newPts[i]));
            lines.push_back(std::move(newLns[i]));
        }
    }


};
