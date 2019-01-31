#include "poly.h"
#include <iomanip>
#include <map>
#include <set>
#include <unordered_set>

#include <GL/glew.h>
#include <GL/gl.h>
#include <GLFW/glfw3.h>
#include <algorithm>
#include <iostream>
#include<string>
#include <cmath>
#include <limits>
#include <cmath>
#include <utility>
#include <iomanip>
#include "linmath.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "intersect.h"
#include <fstream>
#include "surfaceintersect.h"

typedef PTR<Object<Parameter>> Real;
typedef PTR<Object<PPoly3>> Poly3;

// The univariate polynomial resulting from substituting for two
// variables in a trivariate polynomial
class Sub2Poly : public Object<PPoly<Parameter>> {
  Poly3 ppoly3;
  int ixyz;
  Real x, y, z;
  PPoly<Parameter> calculate ();
public:
  // ixyz = 0 or 1 or 2 means substituting for y,z or x,z or x,y
  Sub2Poly (Poly3 ppoly3, int ixyz, Real xyz[3])
    : ppoly3(ppoly3), ixyz(ixyz), x(xyz[0]), y(xyz[1]), z(xyz[2]) {}
};

PPoly<Parameter> Sub2Poly::calculate () {
    return ppoly3->get().substitute2(ixyz, x->get(), y->get(), z->get());
  }

// a + (i/n) (b - a)
class IofN : public Object<Parameter> {
  int i, n;
  Real a, b;
public:
  IofN (int i, int n, Real a, Real b)
    : i(i), n(n), a(a), b(b) {}
  Parameter calculate () { return a->get() + i * ((b->get() - a->get()) / n); }
};

Roots edgeRoots (Poly3 poly3, Real xyz[3], int istep, Real step) {
  ObjPTR<PPoly<Parameter>> poly = new Sub2Poly(poly3, istep, xyz);
  return PolySolver(poly).getRoots(xyz[istep], step);
}


class Square {
public:
    int xyz[3];
    int istep;
    Square() {

    }
    Square(int x, int y, int z, int istep) {
        xyz[0] = x;
        xyz[1] = y;
        xyz[2] = z;
        this->istep = istep;
    }

};

bool operator==(const Square& lhs, const Square& rhs) {
    return lhs.xyz[0] == rhs.xyz[0] &&
           lhs.xyz[1] == rhs.xyz[1] &&
           lhs.xyz[2] == rhs.xyz[2] &&
           lhs.istep == rhs.istep;
}


class GridRoots {
public:
  int nSteps[3];
  Real minXYZ[3];


  void setMinXYZ (Real xyz[3]) {
    for (int i = 0; i < 3; i++)
      minXYZ[i] = xyz[i];
  }

  Real maxXYZ[3];

  void setMaxXYZ (Real xyz[3]) {
    for (int i = 0; i < 3; i++)
      maxXYZ[i] = xyz[i];
  }

  vector<Real> allXYZ[3];

  void setSteps (int nSteps[3]) {
    for (int i = 0; i < 3; i++) {
      this->nSteps[i] = nSteps[i];
      allXYZ[i].push_back(minXYZ[i]);
      for (int j = 1; j < nSteps[i]; j++)
        allXYZ[i].push_back(new IofN(j, nSteps[i], minXYZ[i], maxXYZ[i]));
      allXYZ[i].push_back(maxXYZ[i]);
    }
  }

  class EdgeRoots {
  public:
    int xyz[3]; // grid "point" such as 7, 3, 4
    // actual grid point (allXYZ[0][xyz[0]], allXYZ[1][xyz[1]], allXYZ[2][xyz[2]])
    int istep; // istep = 0, 1, 2 means edge in x, y, z direction
    Roots roots; // intersections on that edge

    EdgeRoots(int xyz[3], int istep, Roots roots)
      : istep(istep), roots(roots) {
      for (int i = 0; i < 3; i++)
        this->xyz[i] = xyz[i];
    }


  };

  vector<EdgeRoots> allRoots;

  void calcEdgeRoots (Poly3 poly3) {
      allRoots.clear();
    int xyz[3], istep;
    Real rxyz[3];
    for (xyz[0] = 0; xyz[0] < allXYZ[0].size(); xyz[0]++) {
      rxyz[0] = allXYZ[0][xyz[0]];
      for (xyz[1] = 0; xyz[1] < allXYZ[1].size(); xyz[1]++) {
        rxyz[1] = allXYZ[1][xyz[1]];
        for (xyz[2] = 0; xyz[2] < allXYZ[2].size(); xyz[2]++) {
          rxyz[2] = allXYZ[2][xyz[2]];
          for (int istep = 0; istep < 3; istep++)
            if (xyz[istep] < allXYZ[istep].size()-1) {
              Roots roots =
                edgeRoots(poly3, rxyz, istep, allXYZ[istep][xyz[istep]+1]);
              if (roots.size() > 0)
                allRoots.push_back(EdgeRoots(xyz, istep, roots));
            }
        }
      }
    }
  }
};
void makeSquare(const GridRoots::EdgeRoots& r, Square s[4]) {
    if (r.istep == 0) {
        s[0] = Square(r.xyz[0], r.xyz[1], r.xyz[2], 1);
        s[1] = Square(r.xyz[0], r.xyz[1], r.xyz[2], 2);
        s[2] = Square(r.xyz[0], r.xyz[1] - 1, r.xyz[2], 1);
        s[3] = Square(r.xyz[0], r.xyz[1], r.xyz[2] - 1, 2);
    }
    else if (r.istep == 1) {
        s[0] = Square(r.xyz[0], r.xyz[1], r.xyz[2], 0);
        s[1] = Square(r.xyz[0], r.xyz[1], r.xyz[2], 2);
        s[2] = Square(r.xyz[0] - 1, r.xyz[1], r.xyz[2], 0);
        s[3] = Square(r.xyz[0], r.xyz[1], r.xyz[2] - 1, 2);
    }
    else {
        s[0] = Square(r.xyz[0], r.xyz[1], r.xyz[2], 0);
        s[1] = Square(r.xyz[0], r.xyz[1], r.xyz[2], 1);
        s[2] = Square(r.xyz[0] - 1, r.xyz[1], r.xyz[2], 0);
        s[3] = Square(r.xyz[0], r.xyz[1]-  1, r.xyz[2], 1);
    }
}

bool neighbors(const GridRoots::EdgeRoots& a, const GridRoots::EdgeRoots& b) {
    Square sa[4];
    Square sb[4];
    makeSquare(a, sa);
    makeSquare(b, sb);

    return sa[0] == sb[0] || sa[0] == sb[1] || sa[0] == sb[2] || sa[0] == sb[3] ||
           sa[1] == sb[0] || sa[1] == sb[1] || sa[1] == sb[2] || sa[1] == sb[3] ||
           sa[2] == sb[0] || sa[2] == sb[1] || sa[2] == sb[2] || sa[2] == sb[3] ||
           sa[3] == sb[0] || sa[3] == sb[1] || sa[3] == sb[2] || sa[3] == sb[3];
}
/*class vertex {
public:
    float x, y, z;
};*/
class Cube {
public:

    bool operator< (const Cube& rhs) const {
      if (this->x < rhs.x) {
        return true;
      }
      else if (this->x > rhs.x) {
        return false;
      }
      else {
        if (this->y < rhs.y) {
          return true;
        }
        else if (this->y > rhs.y) {
          return false;
        }
        else {
          return this->z < rhs.z;
        }
      }
    }

    int x, y, z;
    Cube(int x1, int y1, int z1):
            x{x1}, y{y1}, z{z1} {}
};



static void error_callback(int error, const char* description)
{
    fprintf(stderr, "Error: %s\n", description);
}

mat4x4 model_mat, proj_mat, mvp;

float deg_x = 0.1;
float deg_y = 0.1;
float deg_z = 0.1;

float scale_factor = 1.0;

float inc_deg = 0.1;
static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
static void mouse(GLFWwindow* window, int button, int action, int mods);
static void motion(GLFWwindow* window, double x, double y);
static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);



class DrawVerticesInCube {

public:
    GridRoots grid;
    map<Cube, vector<GridRoots::EdgeRoots>> verticesInCube;
    map<Cube, GLuint> vboForCube;
    GLuint vertex_shader, fragment_shader, program;
    GLint mvp_location;
    GLint vpos_location;
    GLint vcol_location;
    GLFWwindow* window;


    glm::mat4 transform;
    Poly3 poly3;
    std::vector<vertex> vertex_all;
    std::vector<vertex> vertex_all_intersect_points;

    map<Poly3,std::vector<vertex>> surfaceVerteces;
    map<Poly1*,std::vector<vertex>> surfaceVerteces_poly1;
    const char* vertex_shader_text =
            "uniform mat4 MVP;\n"
            "attribute vec3 vCol;\n"
            "attribute vec3 vPos;\n"
            "varying vec3 color;\n" 
            "void main()\n"
            "{\n"
            "    gl_Position = MVP * vec4(vPos, 1.0);\n"
            "    color = vCol;//vec3(1.f,0.f,0.f);\n"
            "}\n";

    const char* fragment_shader_text =
            "varying vec3 color;\n"
            "void main()\n"
            "{\n"
            "    gl_FragColor = vec4(color, 1.0);\n"
            "}\n";

    DrawVerticesInCube() = default;


    void setGrid(const GridRoots& grid) {
        this->grid = grid;
//        grid.setMaxXYZ()
    }

    void setMaxXYZ(Real xyz[3]) {
        grid.setMaxXYZ(xyz);
    }
    void setMinXYZ(Real xyz[3]) {
        grid.setMinXYZ(xyz);
    }

    void setSteps(int nSteps[3]) {
        this->grid.setSteps(nSteps);
    }

    void calcEdgeRoots() {
        this->grid.calcEdgeRoots(poly3);
    }

    void updateVerticesInCube() {
        verticesInCube = findEdgeRootsInCube();
    }

    void setPoly3(Poly3 poly3) {
        this->poly3 = poly3;
    }

    void checkIntersectionPointsHide()
    {
        std::vector<vertex> hidder_verteces;
        zbuffer(transform,vertex_all,vertex_all_intersect_points,hidder_verteces);
        printf("%d points hidden\n",(int)hidder_verteces.size());
    }

    void AddIntersectionTriangles()
    {
        int num=surfaceVerteces.size();
        map<Poly3,std::vector<vertex>>::iterator iter=surfaceVerteces.begin();
        map<Poly3,std::vector<vertex>>::iterator iter1;
        while(iter!=surfaceVerteces.end())
        {
            iter1 = iter;
            ++iter1;
            if(iter1!=surfaceVerteces.end())
            {
                //calc two surface intersections
                auto& verteces1 = *iter;
                auto& verteces2 = *iter1;
                std::vector<vertex> intersectTri;
                std::vector<vertex> intersectPoints;
                SurfaceIntersection(verteces1.second,verteces2.second,intersectTri,intersectPoints);

                int num=intersectTri.size();
                for(int i=0;i<num;++i)
                {
                    vertex color;
                    color.x=1;color.y=0.5;color.z=0.5;
                    vertex_all.push_back(color);
                    vertex_all.push_back(intersectTri[i]);
                }
                vertex_all_intersect_points.insert(vertex_all_intersect_points.end(),
                                                   intersectPoints.begin(),intersectPoints.end());
                intersectTri.clear();
                intersectPoints.clear();
            }
            ++iter;
        }

    }

   void AddIntersectionTriangles2()
    {
        int num=surfaceVerteces_poly1.size();
        map<Poly1*,std::vector<vertex>>::iterator iter=surfaceVerteces_poly1.begin();
        map<Poly1*,std::vector<vertex>>::iterator iter1;
        while(iter!=surfaceVerteces_poly1.end())
        { 
            iter1 = iter;
            ++iter1;
            while(iter1!=surfaceVerteces_poly1.end())
            {
                //calc two surface intersections
                auto& verteces1 = *iter;
                auto& verteces2 = *iter1;
                std::vector<vertex> intersectTri;
                std::vector<vertex> intersectPoints;
                SurfaceIntersection(verteces1.second,verteces2.second,intersectTri,intersectPoints);
                int num=intersectTri.size();
                for(int i=0;i<num;++i)
                {
                    vertex color;
                    color.x=1;color.y=0.5;color.z=0.5;
                    vertex_all.push_back(color);
                    vertex_all.push_back(intersectTri[i]);
                }
				num = intersectPoints.size();
				for (int i = 0; i<num; ++i)
				{
					vertex color;
					color.x = 1.0; color.y = 0.0; color.z = 0.0;
					vertex_all_intersect_points.push_back(color);
					vertex_all_intersect_points.push_back(intersectPoints[i]);
				}
                //vertex_all_intersect_points.insert(vertex_all_intersect_points.end(),
                //                                   intersectPoints.begin(),intersectPoints.end());
                intersectTri.clear();
                intersectPoints.clear(); 
                ++iter1;
            }
            ++iter;
        }

    }

     void AddSurface2(const Poly1& poly1)
     {
         std::vector<Triangle> vecTriangles;
         std::vector<Triangle> vecColors;
         CalcIntersectTriangles(poly1,vecTriangles,vecColors);
         //ComputeIntersectTriangles(vecTriangles,vecColors);

	std::vector<vertex> vecPolyVertex;
         for(int i=0;i<vecTriangles.size();++i)
         {            
             vertex_all.push_back(vecColors[i].v1);
             vertex_all.push_back(vecTriangles[i].v1);
              
             vertex_all.push_back(vecColors[i].v2);
             vertex_all.push_back(vecTriangles[i].v2);
             
             vertex_all.push_back(vecColors[i].v3);
             vertex_all.push_back(vecTriangles[i].v3);

	     vecPolyVertex.push_back(vecTriangles[i].v1);
	     vecPolyVertex.push_back(vecTriangles[i].v2);
	     vecPolyVertex.push_back(vecTriangles[i].v3);	     
         }
	 
	 surfaceVerteces_poly1[(Poly1*)&poly1] = vecPolyVertex; 
	 
         return;
     }

    void AddSurface(const GridRoots& grid, const Poly3& poly3)
    {

        this->grid=grid;
        this->poly3=poly3;
        verticesInCube.clear();

        this->grid.calcEdgeRoots(poly3);
        verticesInCube = findEdgeRootsInCube();
        updateVerticesInCube();
         PPoly3 ppoly3 = poly3->getCurrentValue();
        for (auto& cube_vertices: verticesInCube)
        {
            auto& cube = cube_vertices.first;
            vector<GridRoots::EdgeRoots> eroots = cube_vertices.second;
            int num=eroots.size();
            int nTri = num-3+1;

            for(int i=0;i<nTri;++i)
               for(int m=i;m<i+3;m++)
            //for (auto r: eroots)
            {
                vertex color;
                vertex v;
                auto r=eroots[m%num];

                v.x = grid.allXYZ[0][r.xyz[0]]->getApprox().mid();
                v.y = grid.allXYZ[1][r.xyz[1]]->getApprox().mid();
                v.z = grid.allXYZ[2][r.xyz[2]]->getApprox().mid();

                for (int k = 0; k <  r.roots.size();++k) {
                    if (r.istep == 0) {
                        v.x = r.roots[k]->getApprox().mid();
                    }
                    else if (r.istep == 1) {
                        v.y = r.roots[k]->getApprox().mid();
                    }
                    else {
                        v.z = r.roots[k]->getApprox().mid();
                    }

                   }
                //printf("%d,%d,%d: \n",r.xyz[0],r.xyz[1],r.xyz[2]);
                //printf("%f,%f,%f\n",v.x,v.y,v.z);

                //smooth shader
                //calc gradient
                PV3 grad = ppoly3.gradient(v.x,v.y,v.z);
                PV3 grad_unit = grad.unit();
                color.x=grad_unit.x.mid()*0.5+0.5;
                color.y=grad_unit.y.mid()*0.5+0.5;
                color.z=grad_unit.z.mid()*0.5+0.5;
                vertex_all.push_back(color);
                vertex_all.push_back(v);

                surfaceVerteces[poly3].push_back(v);
            }
        }

    }

    explicit DrawVerticesInCube(const GridRoots& grid, const Poly3& poly3) {
        this->grid = grid;

        verticesInCube = findEdgeRootsInCube();


        initGLFW();

        bufferData();
        initShader();


        runGLFW();
    }


    void runGLFW() {
        while (!glfwWindowShouldClose(window))
        {
            float ratio;
            int width, height;

            glfwGetFramebufferSize(window, &width, &height);
            ratio = width / (float) height;

            glViewport(0, 0, width, height);

            glClearColor(0.0, 0.0, 0.0, 1.0);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            //glUseProgram(program);
            //glUniformMatrix4fv(mvp_location, 1, GL_FALSE, (const GLfloat*) mvp);

            //glPolygonMode(GL_FRONT, GL_LINE);
            //glPolygonMode(GL_BACK, GL_LINE);
            draw2();
			draw3();
            //glPolygonMode(GL_FRONT, GL_FILL);
            //glPolygonMode(GL_BACK, GL_FILL);

            glfwSwapBuffers(window);
            glfwPollEvents();
        }

        glfwDestroyWindow(window);

        glfwTerminate();
    }
    void initGLFW() {
        if (!glfwInit())
            exit(EXIT_FAILURE);
        glfwSetErrorCallback(error_callback);

        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);

        window = glfwCreateWindow(640, 480, "Simple example", NULL, NULL);
        glfwSetMouseButtonCallback(window, mouse);	// callback for mouse click inputs
        glfwSetCursorPosCallback(window, motion);		// callback for mouse movements
        glfwSetScrollCallback(window, scroll_callback);//scole call back

        if (!window)
        {
            glfwTerminate();
            exit(EXIT_FAILURE);
        }

        glfwSetKeyCallback(window, key_callback);
      //  glfwSetScrollCallback(window, scroll_callback);


        glfwMakeContextCurrent(window);
        if (glewInit() != GLEW_OK) {
            fprintf(stderr, "Failed to initialize GLEW\n");
            exit(EXIT_FAILURE);
        }

        glfwSwapInterval(1);


        float ratio;
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        ratio = width / (float) height;

        mat4x4_identity(model_mat);
        mat4x4_rotate_X(model_mat, model_mat, deg_x);
        mat4x4_rotate_Y(model_mat, model_mat, deg_y);
        mat4x4_rotate_Z(model_mat, model_mat, deg_z);
//        mat4x4_ortho(proj_mat, -ratio, ratio, -1.0f, 1.0f, -1.0f, 1.0f);
//        mat4x4_ortho(proj_mat, -1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f);
        mat4x4_ortho(proj_mat, -scale_factor, scale_factor, -scale_factor, scale_factor, -scale_factor, scale_factor);
        mat4x4_mul(mvp, proj_mat, model_mat);


		glDisable(GL_CULL_FACE);
		glEnable(GL_DEPTH_TEST);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDepthFunc(GL_LEQUAL);

    }

    void initShader() {
        vertex_shader = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertex_shader, 1, &vertex_shader_text, NULL);
        glCompileShader(vertex_shader);

        fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragment_shader, 1, &fragment_shader_text, NULL);
        glCompileShader(fragment_shader);

        program = glCreateProgram();
        glAttachShader(program, vertex_shader);
        glAttachShader(program, fragment_shader);
        glLinkProgram(program);

        mvp_location = glGetUniformLocation(program, "MVP");
        vpos_location = glGetAttribLocation(program, "vPos");
        vcol_location = glGetAttribLocation(program, "vCol");

        transform = glm::mat4(1.0f);
        glUseProgram(program);
        glUniformMatrix4fv(mvp_location, 1, GL_FALSE, glm::value_ptr(transform));
    }

    void draw() {
//        glEnableVertexAttribArray(vpos_location);
//        glVertexAttribPointer(vpos_location, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), nullptr);
//        glDrawArrays(GL_POINTS, 1, all_vertices.size());

        for (auto& cube_vertices: verticesInCube) {
             auto& cube = cube_vertices.first;
            auto& vertices = cube_vertices.second;    
            auto& vbo = vboForCube[cube];

            size_t offset = sizeof(vertex);
            glBindBuffer(GL_ARRAY_BUFFER, vbo);
            glEnableVertexAttribArray(vpos_location);
            glVertexAttribPointer(vpos_location, 3, GL_FLOAT, GL_FALSE, 2*sizeof(vertex), (GLvoid*)offset);

            glBindBuffer(GL_ARRAY_BUFFER, vbo);
            glEnableVertexAttribArray(vcol_location);
            glVertexAttribPointer(vcol_location, 3, GL_FLOAT, GL_FALSE, 2*sizeof(vertex), nullptr);

            if (vertices.size() >= 4) {
               // glDrawArrays(GL_LINE_STRIP, 0, vertices.size());

            }
            glDrawArrays(GL_TRIANGLES, 0, 3);
            if (vertices.size() == 3) {
                //glDrawArrays(GL_TRIANGLES, 0, 3);
            }
            else if (vertices.size() == 2) {
//                glDrawArrays(GL_LINE, 0, 2);
            }
            else if (vertices.size() == 1) {
//                glDrawArrays(GL_POINT, 0, 1);
            }

        }
    }

    GLuint all_buffer;
    vector<vertex> all_vertices;
    GLuint abuffer;
    void bufferData() {
        //glGenVertexArrays(1, &abuffer);
        //glBindVertexArray(abuffer);
        PPoly3 ppoly3 = poly3->getCurrentValue();
         for (auto& cube_vertices: verticesInCube)
        {
             auto& cube = cube_vertices.first;
            vector<GridRoots::EdgeRoots> eroots = cube_vertices.second;
            vertex color;
            vertex v;
            vector<vertex> vertices;

            int num=eroots.size();
            int nTri = num-3+1;
            for(int i=0;i<nTri;++i)
                for(int m=i;m<i+3;++m)
            //for (auto r: eroots)
            {
                auto r=eroots[m];
                v.x = grid.allXYZ[0][r.xyz[0]]->getApprox().mid();
                v.y = grid.allXYZ[1][r.xyz[1]]->getApprox().mid();
                v.z = grid.allXYZ[2][r.xyz[2]]->getApprox().mid();

                for (int k = 0; k <  r.roots.size();++k) {
                    if (r.istep == 0) {
                        v.x = r.roots[k]->getApprox().mid();
                    }
                    else if (r.istep == 1) {
                        v.y = r.roots[k]->getApprox().mid();
                    }
                    else {
                        v.z = r.roots[k]->getApprox().mid();
                    }

                   }
                //printf("%d,%d,%d: ",r.xyz[0],r.xyz[1],r.xyz[2]);
                //printf("%f,%f,%f\n",v.x,v.y,v.z);

                //smooth shader
                //calc gradient
                PV3 grad = ppoly3.gradient(v.x,v.y,v.z);
                PV3 grad_unit = grad.unit();
                color.x=grad_unit.x.mid();color.y=grad_unit.y.mid();color.z=grad_unit.z.mid();

                vertices.push_back(color);
                vertices.push_back(v);
            }

            GLuint buffer;
            glGenBuffers(1, &buffer);
            glBindBuffer(GL_ARRAY_BUFFER, buffer);
            glBufferData(GL_ARRAY_BUFFER, sizeof(vertex)*vertices.size(), vertices.data(), GL_STATIC_DRAW);

//            glEnableVertexAttribArray(vpos_location);
//            glVertexAttribPointer(vpos_location, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), nullptr);
//
            vboForCube[cube] = buffer;
        }


//        cout << sizeof(vertex) << endl;
//        cout << all_vertices.size() << endl;
//        glGenBuffers(1, &all_buffer);
//        glBindBuffer(GL_ARRAY_BUFFER, all_buffer);
//        glBufferData(GL_ARRAY_BUFFER, sizeof(vertex)*all_vertices.size(), all_vertices.data(), GL_STATIC_DRAW);

    }
    GLuint vbo;
    void draw2() {
        glUseProgram(program);
        glUniformMatrix4fv(mvp_location, 1, GL_FALSE, glm::value_ptr(transform));

        size_t offset = sizeof(vertex);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glEnableVertexAttribArray(vpos_location);
        glVertexAttribPointer(vpos_location, 3, GL_FLOAT, GL_FALSE, 2*sizeof(vertex), (GLvoid*)offset);

        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glEnableVertexAttribArray(vcol_location);
        glVertexAttribPointer(vcol_location, 3, GL_FLOAT, GL_FALSE, 2*sizeof(vertex), nullptr);
        glDrawArrays( GL_TRIANGLES, 0, vertex_all.size()/2);
		glUseProgram(0);
    }

    void bufferData2() {

        GLuint buffer;
        glGenBuffers(1, &buffer);
        glBindBuffer(GL_ARRAY_BUFFER, buffer);
        int si=vertex_all.size();
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertex)*vertex_all.size(), &vertex_all[0], GL_STATIC_DRAW);

        vbo=buffer;
    }

	GLuint vbo3;
	void draw3() {
		glUseProgram(program);
		glUniformMatrix4fv(mvp_location, 1, GL_FALSE, glm::value_ptr(transform));

		size_t offset = sizeof(vertex);
		glBindBuffer(GL_ARRAY_BUFFER, vbo3);
		glEnableVertexAttribArray(vpos_location);
		glVertexAttribPointer(vpos_location, 3, GL_FLOAT, GL_FALSE, 2 * sizeof(vertex), (GLvoid*)offset);

		glBindBuffer(GL_ARRAY_BUFFER, vbo3);
		glEnableVertexAttribArray(vcol_location);
		glVertexAttribPointer(vcol_location, 3, GL_FLOAT, GL_FALSE, 2 * sizeof(vertex), nullptr);
		glDrawArrays(GL_POINTS, 0, vertex_all_intersect_points.size() / 2);
		glUseProgram(0);
	}

	void bufferData3() {

		GLuint buffer;
		glGenBuffers(1, &buffer);
		glBindBuffer(GL_ARRAY_BUFFER, buffer);
		int si = vertex_all_intersect_points.size();
		glBufferData(GL_ARRAY_BUFFER, sizeof(vertex)*vertex_all_intersect_points.size(), &vertex_all_intersect_points[0], GL_STATIC_DRAW);

		vbo3 = buffer;
	}

    map<Cube, vector<GridRoots::EdgeRoots>> findEdgeRootsInCube() {
        map<Cube, vector<GridRoots::EdgeRoots>> edgerootsInCube;
        for (int i = 0; i < grid.allRoots.size(); i++) {
            auto& r = grid.allRoots[i];

            if (r.istep == 0) {
                for (int k = 0; k < r.roots.size(); ++k) {
                    if (r.xyz[0] < grid.nSteps[0] &&
                        r.xyz[1] < grid.nSteps[1] &&
                        r.xyz[2] < grid.nSteps[2]) {
                        auto c1 = Cube(r.xyz[0], r.xyz[1], r.xyz[2]);
                        edgerootsInCube[c1].push_back(r);
                    }

                    if (r.xyz[0] < grid.nSteps[0] &&
                        r.xyz[1] > 0 &&
                        r.xyz[2] < grid.nSteps[2]) {
                        auto c2 = Cube(r.xyz[0], r.xyz[1]-1, r.xyz[2]);
                        edgerootsInCube[c2].push_back(r);
                    }

                    if (r.xyz[0] < grid.nSteps[0] &&
                        r.xyz[1] < grid.nSteps[1] &&
                        r.xyz[2] > 0) {
                        auto c3 = Cube(r.xyz[0], r.xyz[1], r.xyz[2]-1);
                        edgerootsInCube[c3].push_back(r);
                    }

                    if (r.xyz[0] < grid.nSteps[0] &&
                        r.xyz[1] > 0 &&
                        r.xyz[2] > 0) {
                        auto c4 = Cube(r.xyz[0], r.xyz[1]-1, r.xyz[2]-1);
                        edgerootsInCube[c4].push_back(r);
                    }
                }
            }
            else if (r.istep == 1) {
                for (int k = 0; k < r.roots.size(); ++k) {
                    if (r.xyz[0] < grid.nSteps[0] &&
                        r.xyz[1] < grid.nSteps[1] &&
                        r.xyz[2] < grid.nSteps[2]) {
                        auto c1 = Cube(r.xyz[0], r.xyz[1], r.xyz[2]);
                        edgerootsInCube[c1].push_back(r);
                    }

                    if (r.xyz[0] > 0 &&
                        r.xyz[1] < grid.nSteps[1] &&
                        r.xyz[2] < grid.nSteps[2]) {
                        auto c2 = Cube(r.xyz[0]-1, r.xyz[1], r.xyz[2]);
                        edgerootsInCube[c2].push_back(r);
                    }

                    if (r.xyz[0] < grid.nSteps[0] &&
                        r.xyz[1] < grid.nSteps[1] &&
                        r.xyz[2] > 0) {
                        auto c3 = Cube(r.xyz[0], r.xyz[1], r.xyz[2]-1);
                        edgerootsInCube[c3].push_back(r);
                    }
                    if (r.xyz[0] > 0 &&
                        r.xyz[1] < grid.nSteps[1] &&
                        r.xyz[2] > 0) {
                        auto c4 = Cube(r.xyz[0]-1, r.xyz[1], r.xyz[2]-1);
                        edgerootsInCube[c4].push_back(r);
                    }
                }
            }
            else {
                for (int k = 0; k < r.roots.size(); ++k) {
                    if (r.xyz[0] < grid.nSteps[0] &&
                        r.xyz[1] < grid.nSteps[1] &&
                        r.xyz[2] < grid.nSteps[2]) {
                        auto c1 = Cube(r.xyz[0], r.xyz[1], r.xyz[2]);
                        edgerootsInCube[c1].push_back(r);
                    }

                    if (r.xyz[0] > 0 &&
                        r.xyz[1] < grid.nSteps[1] &&
                        r.xyz[2] < grid.nSteps[2]) {
                        auto c2 = Cube(r.xyz[0]-1, r.xyz[1], r.xyz[2]);
                        edgerootsInCube[c2].push_back(r);
                    }

                    if (r.xyz[0] < grid.nSteps[0] &&
                        r.xyz[1] > 0 &&
                        r.xyz[2] < grid.nSteps[2]) {
                        auto c3 = Cube(r.xyz[0], r.xyz[1]-1, r.xyz[2]);
                        edgerootsInCube[c3].push_back(r);
                    }
                    if (r.xyz[0] > 0 && r.xyz[1] > 0
                        && r.xyz[2] < grid.nSteps[2]) {
                        auto c4 = Cube(r.xyz[0]-1, r.xyz[1]-1, r.xyz[2]);
                        edgerootsInCube[c4].push_back(r);
                    }
                }
            }
        }


        for (auto cube_edge: edgerootsInCube) {
            auto& edges= cube_edge.second;

            for (int i = 0; i < edges.size(); ++i) {
                for (int j = i+1; j < edges.size(); ++j) {
                    if (neighbors(edges[i], edges[j])) {
                        auto temp = GridRoots::EdgeRoots(edges[j].xyz, edges[j].istep, edges[j].roots);
                        edges[j] = edges[i+1];
                        edges[i+1] = temp;
                        break;
                    }
                }
            }

//            cout << cube_edge.second.size() << endl;
        }
        return edgerootsInCube;
    };
};


auto drawer = DrawVerticesInCube();
// Set number of steps to 10 for each dimension.
int nSteps[3] = { 10, 10, 10 };
// Set lower bounds to (-0.5, -0.5, -0.5)
Real lo = new InputObject<Parameter>(Parameter::constant(-1.0));
Real los[3] = { lo, lo, lo };


// Set upper bounds to (0.5, 0.5, 0.5)
Real hi = new InputObject<Parameter>(Parameter::constant(1.0));
Real his[3] = { hi, hi, hi };


double posx;
double posy;
int button_mode=-1;

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    float ratio;
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);
    ratio = width / (float) height;

    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    else if (key == GLFW_KEY_LEFT && action == GLFW_PRESS) {
        mat4x4_identity(model_mat);
        mat4x4_rotate_X(model_mat, model_mat, deg_x);
        mat4x4_rotate_Y(model_mat, model_mat, deg_y);
//            mat4x4_rotate_Z(model_mat, model_mat, 20);
//        mat4x4_ortho(proj_mat, -ratio, ratio, -1.0f, 1.0f, -1.0f, 1.0f);
//        mat4x4_ortho(proj_mat, -1.0f, -1.0f, -1.0f, 1.0f, -1.0f, 1.0f);
        mat4x4_ortho(proj_mat, -scale_factor, scale_factor, -scale_factor, scale_factor, -scale_factor, scale_factor);
        mat4x4_mul(mvp, proj_mat, model_mat);
        deg_y += inc_deg;
    }
    else if (key == GLFW_KEY_RIGHT && action == GLFW_PRESS) {
        mat4x4_identity(model_mat);
        mat4x4_rotate_X(model_mat, model_mat, deg_x);
        mat4x4_rotate_Y(model_mat, model_mat, -deg_y);
//            mat4x4_rotate_Z(model_mat, model_mat, 20);
//        mat4x4_ortho(proj_mat, -ratio, ratio, -1.0f, 1.0f, -1.0f, 1.0f);
//        mat4x4_ortho(proj_mat, -1.0f, -1.0f, -1.0f, 1.0f, -1.0f, 1.0f);
        mat4x4_ortho(proj_mat, -scale_factor, scale_factor, -scale_factor, scale_factor, -scale_factor, scale_factor);
        mat4x4_mul(mvp, proj_mat, model_mat);
        deg_y -= inc_deg;
    }
    else if (key == GLFW_KEY_UP && action == GLFW_PRESS) {
        mat4x4_identity(model_mat);
        mat4x4_rotate_X(model_mat, model_mat, deg_x);
        mat4x4_rotate_Y(model_mat, model_mat, deg_y);
//            mat4x4_rotate_Z(model_mat, model_mat, 20);
//        mat4x4_ortho(proj_mat, -ratio, ratio, -1.0f, 1.0f, -1.0f, 1.0f);
//        mat4x4_ortho(proj_mat, -1.0f, -1.0f, -1.0f, 1.0f, -1.0f, 1.0f);
        mat4x4_ortho(proj_mat, -scale_factor, scale_factor, -scale_factor, scale_factor, -scale_factor, scale_factor);
        mat4x4_mul(mvp, proj_mat, model_mat);
        deg_x += inc_deg;
    }
    else if (key == GLFW_KEY_DOWN && action == GLFW_PRESS) {
        mat4x4_identity(model_mat);
        mat4x4_rotate_X(model_mat, model_mat, -deg_x);
        mat4x4_rotate_Y(model_mat, model_mat, deg_y);
//      mat4x4_rotate_Z(model_mat, model_mat, 20);
//        mat4x4_ortho(proj_mat, -ratio, ratio, -1.0f, 1.0f, -1.0f, 1.0f);
//        mat4x4_ortho(proj_mat, -1.0f, -1.0f, -1.0f, 1.0f, -1.0f, 1.0f);
        mat4x4_ortho(proj_mat, -scale_factor, scale_factor, -scale_factor, scale_factor, -scale_factor, scale_factor);
        mat4x4_mul(mvp, proj_mat, model_mat);
        deg_x -= inc_deg;
    }
    else if (key == GLFW_KEY_U && action == GLFW_PRESS) {
        scale_factor *= 1.1;
        mat4x4_ortho(proj_mat, -scale_factor, scale_factor, -scale_factor, scale_factor, -scale_factor, scale_factor);
        mat4x4_mul(mvp, proj_mat, model_mat);

        mat4x4 modelmat_inv;
        mat4x4_invert(modelmat_inv, model_mat);
        vec4 r1, r2;
        vec4 a, b;
        a[0] = -scale_factor;
        a[1] = -scale_factor;
        a[2] = -scale_factor;
        a[3] = 1;
        b[0] = scale_factor;
        b[1] = scale_factor;
        b[2] = scale_factor;
        b[3] = 1;

        mat4x4_mul_vec4(r1, modelmat_inv, a);
        mat4x4_mul_vec4(r2, modelmat_inv, b);

        cout << scale_factor << ' ' << r1[0] << ' ' << r1[1] << ' ' << r1[2] << ' ' << r1[3] << endl;
        cout << scale_factor << ' ' << r2[0] << ' ' << r2[1] << ' ' << r2[2] << ' ' << r2[3] << endl;

        auto l = min(min(r1[0], r1[1]), r1[2]);
        auto h = max(max(r2[0], r2[1]), r2[2]);
//        mat4x4_ortho(proj_mat, -ratio, ratio, l, h, l, h);
//        mat4x4_mul(mvp, proj_mat, model_mat);
//
//
//        los[0] = new InputObject<Parameter>(Parameter::constant(2*r1[0]));
//        los[1] = new InputObject<Parameter>(Parameter::constant(2*r1[1]));
//        los[2] = new InputObject<Parameter>(Parameter::constant(2*r1[2]));
//
//        his[0] = new InputObject<Parameter>(Parameter::constant(2*r2[0]));
//        his[1] = new InputObject<Parameter>(Parameter::constant(2*r2[1]));
//        his[2] = new InputObject<Parameter>(Parameter::constant(2*r2[2]));
//
//        los[0] = new InputObject<Parameter>(Parameter::constant(2*l));
//        los[1] = new InputObject<Parameter>(Parameter::constant(2*l));
//        los[2] = new InputObject<Parameter>(Parameter::constant(2*l));
//
//        his[0] = new InputObject<Parameter>(Parameter::constant(2*h));
//        his[1] = new InputObject<Parameter>(Parameter::constant(2*h));
//        his[2] = new InputObject<Parameter>(Parameter::constant(2*h));

        los[0] = new InputObject<Parameter>(Parameter::constant(-scale_factor));
        los[1] = new InputObject<Parameter>(Parameter::constant(-scale_factor));
        los[2] = new InputObject<Parameter>(Parameter::constant(-scale_factor));

        his[0] = new InputObject<Parameter>(Parameter::constant(scale_factor));
        his[1] = new InputObject<Parameter>(Parameter::constant(scale_factor));
        his[2] = new InputObject<Parameter>(Parameter::constant(scale_factor));


        GridRoots grid;
        grid.setMinXYZ(los);
        grid.setMaxXYZ(his);

        grid.setSteps(nSteps);
        drawer.setGrid(grid);

        drawer.verticesInCube.clear();
        drawer.vboForCube.clear();

        drawer.calcEdgeRoots();
        drawer.updateVerticesInCube();
//        cout << drawer.grid.allRoots.size() << ' ' << drawer.verticesInCube.size() << endl;
        drawer.bufferData();

    }
    else if (key == GLFW_KEY_D && action == GLFW_PRESS) {
        scale_factor /= 1.1;

        mat4x4_ortho(proj_mat, -scale_factor, scale_factor, -scale_factor, scale_factor, -scale_factor, scale_factor);
        mat4x4_mul(mvp, proj_mat, model_mat);

        mat4x4 modelmat_inv;
        mat4x4_invert(modelmat_inv, model_mat);
        vec4 r1, r2;
        vec4 a, b;
        a[0] = -scale_factor;
        a[1] = -scale_factor;
        a[2] = -scale_factor;
        a[3] = 1;
        b[0] = scale_factor;
        b[1] = scale_factor;
        b[2] = scale_factor;
        b[3] = 1;

        mat4x4_mul_vec4(r1, modelmat_inv, a);
        mat4x4_mul_vec4(r2, modelmat_inv, b);

        cout << scale_factor << ' ' << r1[0] << ' ' << r1[1] << ' ' << r1[2] << ' ' << r1[3] << endl;
        cout << scale_factor << ' ' << r2[0] << ' ' << r2[1] << ' ' << r2[2] << ' ' << r2[3] << endl;

        auto l = min(min(r1[0], r1[1]), r1[2]);
        auto h = max(max(r2[0], r2[1]), r2[2]);
//        mat4x4_ortho(proj_mat, -ratio, ratio, l, h, l, h);
//        mat4x4_mul(mvp, proj_mat, model_mat);

//
//        los[0] = new InputObject<Parameter>(Parameter::constant(2*r1[0]));
//        los[1] = new InputObject<Parameter>(Parameter::constant(2*r1[1]));
//        los[2] = new InputObject<Parameter>(Parameter::constant(2*r1[2]));
//
//        his[0] = new InputObject<Parameter>(Parameter::constant(2*r2[0]));
//        his[1] = new InputObject<Parameter>(Parameter::constant(2*r2[1]));
//        his[2] = new InputObject<Parameter>(Parameter::constant(2*r2[2]));

//
//        los[0] = new InputObject<Parameter>(Parameter::constant(2*l));
//        los[1] = new InputObject<Parameter>(Parameter::constant(2*l));
//        los[2] = new InputObject<Parameter>(Parameter::constant(2*l));
//
//        his[0] = new InputObject<Parameter>(Parameter::constant(2*h));
//        his[1] = new InputObject<Parameter>(Parameter::constant(2*h));
//        his[2] = new InputObject<Parameter>(Parameter::constant(2*h));

        los[0] = new InputObject<Parameter>(Parameter::constant(-scale_factor));
        los[1] = new InputObject<Parameter>(Parameter::constant(-scale_factor));
        los[2] = new InputObject<Parameter>(Parameter::constant(-scale_factor));

        his[0] = new InputObject<Parameter>(Parameter::constant(scale_factor));
        his[1] = new InputObject<Parameter>(Parameter::constant(scale_factor));
        his[2] = new InputObject<Parameter>(Parameter::constant(scale_factor));

        GridRoots grid;
        grid.setMinXYZ(los);
        grid.setMaxXYZ(his);
        grid.setSteps(nSteps);
        drawer.setGrid(grid);

        drawer.verticesInCube.clear();
        drawer.vboForCube.clear();

        drawer.calcEdgeRoots();
        drawer.updateVerticesInCube();
//        cout << drawer.grid.allRoots.size() << ' ' << drawer.verticesInCube.size() << endl;
        drawer.bufferData();

//        cout<< l << ' ' << h << endl;
//        mat4x4_ortho(proj_mat, -ratio, ratio, l, h, l, h);
//        mat4x4_mul(mvp, proj_mat, model_mat);
    }
}


static void rotate(float dx,float dy)
{
   drawer.transform = glm::rotate(drawer.transform , dy, glm::vec3(0.0f,1.0f,0.0f));//旋转
   drawer.transform = glm::rotate(drawer.transform ,dx, glm::vec3(1.0f,0.0f, 0.0f));//旋转
   glUniformMatrix4fv(drawer.mvp_location, 1, GL_FALSE, glm::value_ptr(drawer.transform ));
}

static void scale(float fScale)
{
    drawer.transform = glm::scale(drawer.transform, glm::vec3(fScale, fScale, fScale));//缩放
    glUniformMatrix4fv(drawer.mvp_location, 1, GL_FALSE, glm::value_ptr(drawer.transform));
}

static void move(float dx,float dy)
{
   double dx1 = drawer.transform[0][0]*dx+drawer.transform[0][1]*dy;
    double dy1 = drawer.transform[1][0]*dx+drawer.transform[1][1]*dy;
    drawer.transform = glm::translate(drawer.transform, glm::vec3(dx1,dy1, 0.0f));//移动
    glUniformMatrix4fv(drawer.mvp_location, 1, GL_FALSE, glm::value_ptr(drawer.transform));
}

static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
if(yoffset<0)
    scale(1.1);
else
    scale(0.9);
}

static void mouse(GLFWwindow* window, int button, int action, int mods)
{

    if (action == GLFW_PRESS)
    {
        glfwGetCursorPos(window, &posx, &posy);
        if(button == GLFW_MOUSE_BUTTON_LEFT)
            button_mode = 1;
        else if(button == GLFW_MOUSE_BUTTON_RIGHT)
            button_mode=2;

    }
    else if(action == GLFW_RELEASE)
    {
       button_mode=-1;
    }
}

static void motion(GLFWwindow* window, double x, double y)
{
    int width, height;
    glfwGetFramebufferSize(drawer.window, &width, &height);
    if(button_mode==1)
    {
        double dx=2*(x-posx)/width;
        double dy=-2*(y-posy)/height;
        move(dx,dy);
        posx=x;
        posy=y;
    }
    else if(button_mode==2)
    {
        double dx=-2*(x-posx)/width;
        double dy=-2*(y-posy)/height;
        rotate(dx,dy);
        posx=x;
        posy=y;
    }
}




void ReadPolyText(char* filename,PPoly3* ppoly3)
{
  std::ifstream fin;
  fin.open(filename);
  string str;

  getline(fin,str);
  getline(fin,str);
  getline(fin,str);
  getline(fin,str);
  int num= atoi(str.c_str())-3;

   const char *d = " ";
  for(int i=0;i<num;++i)
  {
       getline(fin,str);
       char *p;
       p = strtok((char*)str.c_str(),d);
       double tmp[4];
       int m=0;
       while(p)
       {
         tmp[m++] = atof(p);
         p=strtok(NULL,d);
       }
       ppoly3->add(tmp[0], tmp[1], tmp[2], Parameter::input(tmp[3]));
  }
  return;
}


int main1 (int argc, char *argv[]) {
  Parameter::enable();

  // The grid and (eventually) its intersections.
  GridRoots grid;

  grid.setMinXYZ(los);
  grid.setMaxXYZ(his);
  grid.setSteps(nSteps);

  // Create a trivariate f(x,y,z)
  // x^4 + y^3 - z - 0.25
  //test: x2+y2+z2=0.5
  PPoly3 ppoly3;
  ppoly3.add(4, 0, 0, Parameter::input(1));
  ppoly3.add(0, 3, 0, Parameter::input(1));
  ppoly3.add(0, 0, 1, Parameter::input(-1));
  ppoly3.add(0, 0, 0, Parameter::input(-0.25));

  Poly3 poly3 = new InputObject<PPoly3>(ppoly3);


  //2 poly
  PPoly3 ppoly3_1;   
  ppoly3_1.add(0, 0, 0, Parameter::input(1.50));
  ppoly3_1.add(2, 0, 0, Parameter::input(-1));
  ppoly3_1.add(0, 1, 0, Parameter::input(-1));
  ppoly3_1.add(0, 0, 1, Parameter::input(-1));
  //ppoly3_1.add(2, 3, 1, Parameter::input(-0.000211234));
  //ppoly3_1.add(3, 1, 2, Parameter::input(0.000566199));
  //ppoly3_1.add(1, 2, 3, Parameter::input(0.000680376));
  Poly3 poly3_1 = new InputObject<PPoly3>(ppoly3_1);

  //3 poly
  PPoly3 ppoly3_2;
  ppoly3_2.add(0, 0, 0, Parameter::input(1.50));
  ppoly3_2.add(1, 0, 0, Parameter::input(-1));
  ppoly3_2.add(0, 1, 0, Parameter::input(-1));
  ppoly3_2.add(0, 0, 1, Parameter::input(-1));
  //ppoly3_2.add(2, 3, 1, Parameter::input(-0.000211234));
  //ppoly3_2.add(3, 1, 2, Parameter::input(0.000566199));
  //ppoly3_2.add(1, 2, 3, Parameter::input(0.000680376));
  Poly3 poly3_2 = new InputObject<PPoly3>(ppoly3_2);


  //24 poly
  PPoly3 ppoly3_3;
  ppoly3_3.add(0, 0, 0, Parameter::input(-0.500004));
  ppoly3_3.add(1, 0, 0, Parameter::input(1));
  ppoly3_3.add(0, 1, 0, Parameter::input(1));
  ppoly3_3.add(0, 0, 1, Parameter::input(-1));
  //ppoly3_3.add(2, 3, 1, Parameter::input(-4.52058e-05));
  //ppoly3_3.add(3, 1, 2, Parameter::input(0.000257742));
  //ppoly3_3.add(1, 2, 3, Parameter::input(0.00010794));
  Poly3 poly3_3 = new InputObject<PPoly3>(ppoly3_3);


  // Intersect the surface f(x,y,z)=0 with the grid edges.
 // grid.calcEdgeRoots(poly3);
  /*drawer.setGrid(grid);
  drawer.setPoly3(poly3);
  drawer.calcEdgeRoots();
  drawer.updateVerticesInCube();
  drawer.initGLFW();
  drawer.bufferData();
  drawer.initShader();
   drawer.runGLFW();*/

  drawer.AddSurface(grid,poly3_2);
  //drawer.AddSurface(grid,poly3);
  //drawer.AddIntersectionTriangles();
  drawer.initGLFW();

  drawer.bufferData2();
  drawer.initShader();

  drawer.runGLFW();
}


void ReadPolyText1(const char* filename,Poly1* poly_)
{
  std::ifstream fin;
  fin.open(filename);
  if(!fin)  
    {      
	printf("open parameter file failed!\n");
	return;
    } 
  string str;

  getline(fin,str);
  getline(fin,str);
  getline(fin,str);
  getline(fin,str);
  int num= atoi(str.c_str())-3;

   const char *d = " ";
  for(int i=0;i<num;++i)
  {
       getline(fin,str);
       char *p;
       p = strtok((char*)str.c_str(),d);
       double tmp[4];
       int m=0;
       while(p)
       {
         tmp[m++] = atof(p);
         p=strtok(NULL,d);
       }
       poly_->Add(tmp[0], tmp[1], tmp[2], tmp[3]);
  }
  return;
}

int main (int argc, char *argv[]) {
  // Create a trivariate f(x,y,z)
  // x^4 + y^3 - z - 0.25
  //test: x2+y2+z2=0.5
  Poly1 ppoly3;
  ppoly3.Add(4, 0, 0,1);
  ppoly3.Add(0, 3, 0, 1);
  ppoly3.Add(0, 0, 1, -1);
  ppoly3.Add(0, 0, 0, -0.25);

  //1 poly
  Poly1 ppoly3_0;
  /*ppoly3_0.Add(0, 0, 0, 1.49999);
  ppoly3_0.Add(1, 0, 0, -1);
  ppoly3_0.Add(0, 1, 0,-1);
  ppoly3_0.Add(0, 0, 1, -1);
  ppoly3_0.Add(2, 3, 1, -0.000211234);
  ppoly3_0.Add(3, 1, 2, 0.000566199);
  ppoly3_0.Add(1, 2, 3, 0.000680376);*/
  ReadPolyText1(("../surf1.txt"),&ppoly3_0);

  //2 poly
  Poly1 ppoly3_1;
  /*ppoly3_1.Add(0, 0, 0, -0.500011);
  ppoly3_1.Add(1, 0, 0, -1);
  ppoly3_1.Add(0, 1, 0,1);
  ppoly3_1.Add(0, 0, 1, 1);
  ppoly3_1.Add(2, 3, 1, 0.000823295);
  ppoly3_1.Add(3, 1, 2, -0.000604897);
  ppoly3_1.Add(1, 2, 3, 0.000596881);*/
  ReadPolyText1(("../surf2.txt"),&ppoly3_1);

  //3 poly
  Poly1 ppoly3_2;
  /*ppoly3_2.Add(0, 0, 0, -0.499995);
  ppoly3_2.Add(1, 0, 0, 1);
  ppoly3_2.Add(0, 1, 0, -1);
  ppoly3_2.Add(0, 0, 1, 1);
  ppoly3_2.Add(2, 3, 1, 0.00053646);
  ppoly3_2.Add(3, 1, 2, -0.00044445);
  ppoly3_2.Add(1, 2, 3, -0.000329554);*/
  ReadPolyText1(("../surf3.txt"),&ppoly3_2);


  drawer.AddSurface2(ppoly3);
  drawer.AddSurface2(ppoly3_0);
  drawer.AddSurface2(ppoly3_1);
  //drawer.AddSurface2(ppoly3_2);
  //drawer.AddIntersectionTriangles2();
  drawer.initGLFW();

  drawer.bufferData2();
  drawer.bufferData3();
  drawer.initShader();

  drawer.runGLFW();
}
