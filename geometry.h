#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <ostream>
#include <utility>
#include <vector>

class Line;

struct Point {
    double x, y;
    Point() : x(0), y(0) {}
    Point(double x, double y) : x(x), y(y) {}
    double distance_to(const Point& other) const;
    double distance_to(const Line& line) const;
    Point point_rotation(const Point& center_of_rotation, double angle);
    Point point_reflection(const Point& center);
    Point point_reflection(const Line& axis);
    Point point_scale(const Point& center, double coefficient);
};

class Line {
  public:
    Line(double k, double b)
        : a_line_coeff_(k), b_line_coeff_(-1), c_line_coeff_(b) {}
    Line(const Point& A, double k)
        : a_line_coeff_(k),
          b_line_coeff_(-1),
          c_line_coeff_(A.y - a_line_coeff_ * A.x) {}
    Line(const Point& A, const Point& B);
    Line perpendicular_from_point(const Point& point) const;
    double get_a_coeff() const {
        return a_line_coeff_;
    }
    double get_b_coeff() const {
        return b_line_coeff_;
    }
    double get_c_coeff() const {
        return c_line_coeff_;
    }
    Point get_point(double x) const;

  private:
    const double EPS = 1e-9;
    double a_line_coeff_, b_line_coeff_, c_line_coeff_;
};

std::ostream& operator<<(std::ostream& os, const Line& line);
std::ostream& operator<<(std::ostream& os, const Point& point);

bool operator==(const Line& lhs, const Line& rhs);
bool operator!=(const Line& lhs, const Line& rhs);
bool operator==(const Point& lhs, const Point& rhs);
bool operator!=(const Point& lhs, const Point& rhs);
double det(double a, double b, double c, double d);
Point center_of_segment(const Point& start, const Point& end);
Point line_intersections(const Line& first, const Line& second);
std::pair<double, double> solvint_quadratic_equation(double coeff_a,
                                                     double coeff_b,
                                                     double coeff_c);
Line rotation_line_about_point(const Line& line, const Point& point,
                               double angle);
bool equality(double a, double b);

class Shape {
  public:
    virtual double perimeter() const = 0;
    virtual double area() const = 0;
    bool operator==(const Shape& another) const;
    bool operator!=(const Shape& another) const;
    bool isCongruentTo(const Shape& another);
    bool isSimilarTo(const Shape& another);
    virtual bool containsPoint(const Point& point) = 0;
    virtual void rotate(const Point& center, double angle) = 0;
    virtual void reflect(const Point& center) = 0;
    virtual void reflect(const Line& axis) = 0;
    virtual void scale(const Point& center, double coefficient) = 0;
    virtual ~Shape() {}
};

class Ellipse : public Shape {
  public:
    Ellipse() = default;
    Ellipse(Point f1, Point f2, double a) : F1_(f1), F2_(f2), a_(a / 2) {}
    std::pair<Point, Point> focuses() const;
    std::pair<Line, Line> directrices() const;
    double eccentricity() const;
    Point center() const;

    double perimeter() const override;
    double area() const override;
    void rotate(const Point& center, double angle) override;
    void reflect(const Point& center) override;
    void reflect(const Line& axis) override;
    void scale(const Point& center, double coefficient) override;
    bool containsPoint(const Point& point) override;

  protected:
    Point F1_, F2_;
    double a_;
    double get_small_axis() const;
};

class Polygon : public Shape {
  public:
    template <typename... Points>
    Polygon(const Points... points) {
        builder(points...);
    }

    Polygon(const std::vector<Point>& vertexes);
    size_t verticesCount() const;
    std::vector<Point> getVertices() const;
    virtual bool isConvex();

    double perimeter() const override;
    double area() const override;
    void rotate(const Point& center, double angle) override;
    void reflect(const Point& center) override;
    void reflect(const Line& axis) override;
    void scale(const Point& center, double coefficient) override;
    bool containsPoint(const Point& point) override;
    double oriented_area() const;

  protected:
    std::vector<Point> vertexes_;
    void builder() {}

    void builder(const Point& point) {
        vertexes_.push_back(point);
    }

    template <typename Head, typename... Tail>
    void builder(const Head& head, const Tail&... tail) {
        vertexes_.push_back(head);
        builder(tail...);
    }
};

std::ostream& operator<<(std::ostream& os, const Polygon& pol);

class Rectangle : public Polygon {
  public:
    Rectangle(const Point& v1, const Point& v2, double k);
    Point center();
    std::pair<Line, Line> diagonals();

    bool isConvex() override;
};

class Circle : public Ellipse {
  public:
    Circle(const Point& A, double r) : Ellipse(A, A, r) {}
    Circle(const Point& A, const Point& B, const Point& C);
    double radius() const;

    double perimeter() const override;
    double area() const override;
};

std::ostream& operator<<(std::ostream& os, const Circle& circle);

class Square : public Rectangle {
  public:
    Square(const Point& A, const Point& B) : Rectangle(A, B, 1.0) {}

    Circle circumscribedCircle();
    Circle inscribedCircle();
};

class Triangle : public Polygon {
  public:
    Triangle(const Point& A, const Point& B, const Point& C)
        : Polygon(A, B, C) {}

    Circle circumscribedCircle();
    Circle inscribedCircle();
    Point centroid();
    Point orthocenter();
    Line EulerLine();
    Circle ninePointsCircle();
    bool isConvex() override;

  private:
    Line get_median(const Point& point);
    Line get_bisector(const Point& point);
};
