#define _MAIN_ 
#ifdef _MAIN_

#include "main.h"

// zSpace Core Headers
#include <headers/zApp/include/zObjects.h>
#include <headers/zApp/include/zFnSets.h>
#include <headers/zApp/include/zViewer.h>

#include <vector>
#include <limits>
#include <cmath>
#include <fstream>  // Include for file output

#include <algorithm>
#include <tuple>

using namespace zSpace;

#define zPt2Vec(pt) Alice::vec(pt.x, pt.y, pt.z)
#define drawPt(pt) drawPoint(zPt2Vec(pt));

// Constants for grid size
const int eX = 50;
const int eY = 50;
const int Numb_Origins = 2;  // Set the number of origins

// Clamp a value between min and max
float clamp(float value, float min, float max) {
    return std::max(min, std::min(max, value));
}

// Map a scalar value in the range [vmin, vmax] to RGB color in the 0-1 range
std::tuple<float, float, float> scalarToRGB(float v, float vmin, float vmax) {
    // Normalize v to a 0-1 range
    float ratio = clamp((v - vmin) / (vmax - vmin), 0.0f, 1.0f);

    // Calculate RGB components in 0-1 range
    float r = std::min(2 * ratio, 1.0f);
    float g = std::min(2 * (1 - ratio), 1.0f);
    float b = 1.0f - std::abs(2 * ratio - 1.0f);

    return std::make_tuple(r, g, b);
}

class RadialField
{
public:
    zPoint AttrOrigins[Numb_Origins]; // Array to store multiple origin points
    zPoint MyGrid[eX][eY];            // 2D array of grid points
    float AttrSize;                   // Field size scaling factor

    // Container to store filtered points
    std::vector<zPoint> filteredPoints;

    // Flag to control one-time printing
    bool printOnce = false;

    // Initialize grid points
    void setupGrid()
    {
        for (int i = 0; i < eX; i++)
        {
            for (int j = 0; j < eY; j++)
            {
                MyGrid[i][j] = zPoint(i, j, 0); // Initialize each grid point at (i, j, k)
            }
        }
    }

    // Function to find points using intersection logic
    void findPointsByColorValueIntersection(float targetColorValue, float tolerance = 0.01)
    {
        filteredPoints.clear(); // Clear previous results

        for (int i = 0; i < eX; i++)
        {
            for (int j = 0; j < eY; j++)
            {
                zPoint GridPt = MyGrid[i][j];
                float maxDistance = 0.0f; // Start with zero for intersection logic

                for (int n = 0; n < Numb_Origins; n++)
                {
                    zVector dir = GridPt - AttrOrigins[n];
                    float distance = dir.length();
                    maxDistance = std::max(maxDistance, distance); // Update to the largest distance
                }

                float fieldSize = AttrSize * 10;
                float colorValue = maxDistance / fieldSize;

                // Check if colorValue is within the target range
                if (std::fabs(colorValue - targetColorValue) <= tolerance)
                {
                    filteredPoints.push_back(GridPt);
                }
            }
        }
    }

    // Draw radial field using intersection logic
    void drawRadialFieldIntersection()
    {
        for (int i = 0; i < eX; i++)
        {
            for (int j = 0; j < eY; j++)
            {
                zPoint GridPt = MyGrid[i][j];
                float maxDistance = 0.0f; // Start with zero for intersection logic

                for (int n = 0; n < Numb_Origins; n++)
                {
                    zVector dir = GridPt - AttrOrigins[n];
                    float distance = dir.length();
                    maxDistance = std::max(maxDistance, distance); // Update to the largest distance
                }

                float fieldSize = AttrSize * 10;
                float r, g, b;
                std::tie(r, g, b) = scalarToRGB(maxDistance, 0, fieldSize);

                glColor3f(r, g, b); // Set color based on max distance
                drawPt(GridPt);     // Draw the point
            }
        }
    }

    void exportPointsToCSV(const std::string& filename)
    {
        std::ofstream file(filename);
        if (file.is_open())
        {
            file << "X,Y,Z\n";

            // Corrected to 2D iteration over MyGrid
            for (int i = 0; i < eX; i++)
            {
                for (int j = 0; j < eY; j++)
                {
                    zPoint pt = MyGrid[i][j];

                    // Ensure the point is valid
                    if (pt.x != -999 && pt.y != -999 && pt.z != -999)
                    {
                        file << std::fixed << std::setprecision(6) << pt.x << "," << pt.y << "," << pt.z << "\n";
                    }
                }
            }

            file.close();
            std::cout << "Points exported successfully to " << filename << std::endl;
        }
        else
        {
            std::cerr << "Unable to open file " << filename << std::endl;
        }
    }
};

// Global instance of RadialField
RadialField rdf1;

void setup()
{
    // Initialize origins of the radial field
    rdf1.AttrOrigins[0] = zPoint(22.5, 25, 0);
    rdf1.AttrOrigins[1] = zPoint(27.5, 25, 0);

    rdf1.AttrSize = 1; // Adjust size scaling

    rdf1.setupGrid(); // Initialize the grid points

    // Set target color value and tolerance
    float targetColorValue = 0.5; // Target color value (between 0 and 1)
    float tolerance = 0.01;       // Tolerance around the target color value
    rdf1.findPointsByColorValueIntersection(targetColorValue, tolerance); // Use intersection method
}

void draw()
{
    // Set background color to grey
    backGround(0);

    // Draw a grid of size 50 units
    drawGrid(50);

    // Set point size to 5
    glPointSize(5);

    // Draw the radial field based on distances to multiple origins
    rdf1.drawRadialFieldIntersection(); // Draw using intersection logic

    // Draw filtered points in a distinct color (e.g., blue)
    glColor3f(0, 0, 1);
    for (const auto& pt : rdf1.filteredPoints)
    {
        drawPt(pt);

        // Print point coordinates once
        if (!rdf1.printOnce)
        {
            std::cout << "Point coordinates: (" << pt.x << ", " << pt.y << ", " << pt.z << ")\n";
        }
    }

    // Set the printOnce flag to true after printing coordinates
    rdf1.printOnce = true;

    // Export points to CSV
    rdf1.exportPointsToCSV("C:\\Code\\Output\\Rain.csv");
}

void update(int value)
{
    // Update function, if needed
}

void keyPress(unsigned char k, int xm, int ym)
{

}

void mousePress(int b, int state, int x, int y)
{
    // Mouse press event handling
}

void mouseMotion(int x, int y)
{
    // Mouse motion event handling 
}

#endif // _MAIN_
