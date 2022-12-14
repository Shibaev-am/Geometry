#include "geometry.h"

bool operator==(const Line& lhs, const Line& rhs) {
    const double EPS = 1e-9;
    double detAB = det(lhs.get_a_coeff(), lhs.get_b_coeff(), rhs.get_a_coeff(),
                       rhs.get_b_coeff());
    return fabs(detAB) < EPS &&
           fabs(det(lhs.get_a_coeff(), lhs.get_c_coeff(), rhs.get_a_coeff(),
                    rhs.get_c_coeff())) < EPS &&
           fabs(det(lhs.get_b_coeff(), lhs.get_c_coeff(), rhs.get_b_coeff(),
                    rhs.get_c_coeff())) < EPS;
}

bool operator!=(const Line& lhs, const Line& rhs) {
    return !(lhs == rhs);
}

bool operator==(const Point& lhs, const Point& rhs) {
    return (equality(lhs.x, rhs.x) && equality(lhs.y, rhs.y));
}
bool operator!=(const Point& lhs, const Point& rhs) {
    return !(lhs == rhs);
}

double det(double a, double b, double c, double d) {
    return a * d - b * c;
}

double Point::distance_to(const Point& other) const {
    return sqrt((x - other.x) * (x - other.x) + (y - other.y) * (y - other.y));
}

double Point::distance_to(const Line& line) const {
    double res = fabs((line.get_a_coeff() * x + line.get_b_coeff() * y +
                       line.get_c_coeff())) /
                 std::sqrt(line.get_a_coeff() * line.get_a_coeff() +
                           line.get_b_coeff() * line.get_b_coeff());
    return res;
}

Point center_of_segment(const Point& start, const Point& end) {
    return {(start.x + end.x) / 2, (start.y + end.y) / 2};
}

Point line_intersections(const Line& first, const Line& second) {
    const double EPS = 1e-9;
    Point ans;
    double detAB = det(first.get_a_coeff(), first.get_b_coeff(),
                       second.get_a_coeff(), second.get_b_coeff());
    if (fabs(detAB) >= EPS) {
        ans.x = -det(first.get_c_coeff(), first.get_b_coeff(),
                     second.get_c_coeff(), second.get_b_coeff()) /
                detAB;
        ans.y = -det(first.get_a_coeff(), first.get_c_coeff(),
                     second.get_a_coeff(), second.get_c_coeff()) /
                detAB;
        return ans;

    } else {
        throw -1;
    }
}

Line::Line(const Point& A, const Point& B)
    : a_line_coeff_(B.y - A.y),
      b_line_coeff_(A.x - B.x),
      c_line_coeff_(A.y * B.x - A.x * B.y) {}

Point Line::get_point(double x) const {
    Point ans;
    if (fabs(b_line_coeff_) > EPS) {
        ans.x = x;
        ans.y = (-c_line_coeff_ - a_line_coeff_ * x) / b_line_coeff_;
    } else {
        ans.x = -c_line_coeff_ / a_line_coeff_;
        ans.y = ans.x + x;
    }
    return ans;
}

std::pair<double, double> solvint_quadratic_equation(double coeff_a,
                                                     double coeff_b,
                                                     double coeff_c) {
    double discriminant = coeff_b * coeff_b - 4 * coeff_a * coeff_c;
    const double EPS = 1e-9;
    if (fabs(coeff_a) < EPS) {
        return {-coeff_c / coeff_b, -coeff_c / coeff_b};
    }
    return {(-coeff_b + sqrt(discriminant)) / (2 * coeff_a),
            (-coeff_b - sqrt(discriminant)) / (2 * coeff_a)};
}

std::ostream& operator<<(std::ostream& os, const Line& line) {
    os << line.get_a_coeff() << " * x + " << line.get_b_coeff() << " * y + "
       << line.get_c_coeff() << " = 0";
    return os;
}

std::ostream& operator<<(std::ostream& os, const Point& point) {
    os << std::fixed << "(" << point.x << ", " << point.y << ")";
    return os;
}

Line Line::perpendicular_from_point(const Point& point) const {
    Line ans = *this;
    if (fabs(b_line_coeff_) <= EPS) {
        ans.a_line_coeff_ = 0;
        ans.b_line_coeff_ = 1;
        ans.c_line_coeff_ = -point.y;
    } else if (fabs(a_line_coeff_) <= EPS) {
        ans.a_line_coeff_ = 1;
        ans.b_line_coeff_ = 0;
        ans.c_line_coeff_ = -point.x;
    } else {
        Line temp(point, b_line_coeff_ / a_line_coeff_);
        return temp;
    }
    return ans;
}

Point Point::point_rotation(const Point& center_of_rotation, double angle) {
    double m_sin = sin(angle), m_cos = cos(angle);
    double new_x = (x - center_of_rotation.x) * m_cos -
                   (y - center_of_rotation.y) * m_sin + center_of_rotation.x;
    double new_y = (x - center_of_rotation.x) * m_sin +
                   (y - center_of_rotation.y) * m_cos + center_of_rotation.y;
    Point new_point(new_x, new_y);
    return new_point;
}

Line rotation_line_about_point(const Line& line, const Point& point,
                               double angle) {
    Point first_point = line.get_point(0).point_rotation(point, angle);
    Point second_point = line.get_point(1).point_rotation(point, angle);
    Line ans(first_point, second_point);
    return ans;
}

Point Point::point_reflection(const Point& center) {
    return Point(2 * center.x - x, 2 * center.y - y);
}

Point Point::point_reflection(const Line& axis) {
    Line perpendicular = axis.perpendicular_from_point(*this);
    Point center = line_intersections(axis, perpendicular);
    return point_reflection(center);
}

Point Point::point_scale(const Point& center, double coefficient) {
    Point ans((x - center.x) * coefficient + center.x,
              (y - center.y) * coefficient + center.y);
    return ans;
}

bool equality(double a, double b) {
    const double EPS = 1e-9;
    return fabs(a - b) <= EPS;
}

double Circle::radius() const {
    return a_;
}

double Circle::perimeter() const {
    return M_PI * 2 * a_;
}
double Circle::area() const {
    return M_PI * a_ * a_;
}

Circle::Circle(const Point& A, const Point& B, const Point& C) : Ellipse() {
    Point centerAB = center_of_segment(A, B);
    Point centerBC = center_of_segment(B, C);
    Line lineAB(A, B);
    Line lineBC(B, C);
    Line first_serper = lineAB.perpendicular_from_point(centerAB);
    Line second_serper = lineBC.perpendicular_from_point(centerBC);
    Point center_of_circle = line_intersections(first_serper, second_serper);
    F1_ = F2_ = center_of_circle;
    a_ = center_of_circle.distance_to(A);
}

std::ostream& operator<<(std::ostream& os, const Circle& circle) {
    os << "( O" << circle.center() << ", R = " << circle.radius() << ")";
    return os;
}

std::pair<Point, Point> Ellipse::focuses() const {
    return {F1_, F2_};
}

std::pair<Line, Line> Ellipse::directrices() const {
    Point p_center = center();
    Line big_axis(F1_, F2_);
    Point start_first_dir;
    Point start_second_dir;
    if (!equality(big_axis.get_b_coeff(), 0)) {
        std::pair<double, double> line = {
            -big_axis.get_a_coeff() / big_axis.get_b_coeff(),
            -big_axis.get_c_coeff() / big_axis.get_b_coeff()};
        double coeff_a = 1 + line.first * line.first;
        double coeff_b =
            -(2 * p_center.x + 2 * line.first * (p_center.y - line.second));
        double coeff_c =
            (p_center.x) * (p_center.x) +
            (p_center.y - line.second) * (p_center.y - line.second) -
            (a_ / eccentricity()) * (a_ / eccentricity());
        std::pair<double, double> start_of_derectices_x =
            solvint_quadratic_equation(coeff_a, coeff_b, coeff_c);
        start_first_dir = {
            start_of_derectices_x.first,
            line.first * start_of_derectices_x.first + line.second};
        start_second_dir = {
            start_of_derectices_x.second,
            line.first * start_of_derectices_x.second + line.second};
    } else {
        double delta = a_ / eccentricity();
        start_first_dir = {p_center.x, p_center.y + delta};
        start_second_dir = {p_center.x, p_center.y - delta};
    }

    Line first_dir = big_axis.perpendicular_from_point(start_first_dir);
    Line second_dir = big_axis.perpendicular_from_point(start_second_dir);
    return {first_dir, second_dir};
}

double Ellipse::eccentricity() const {
    return (F1_.distance_to(F2_) / (2 * a_));
}

Point Ellipse::center() const {
    return center_of_segment(F1_, F2_);
}

double Ellipse::perimeter() const {
    double b = get_small_axis();
    double temp1 = 3 * ((a_ - b) / (a_ + b)) * ((a_ - b) / (a_ + b));
    double temp2 = 10 + sqrt(4 - temp1);
    return M_PI * (a_ + b) * (1 + temp1 / temp2);
}

double Ellipse::area() const {
    return M_PI * a_ * get_small_axis();
}

void Ellipse::rotate(const Point& center, double angle) {
    F1_ = F1_.point_rotation(center, angle / 180 * M_PI);
    F2_ = F2_.point_rotation(center, angle / 180 * M_PI);
}

void Ellipse::reflect(const Point& center) {
    F1_ = F1_.point_reflection(center);
    F2_ = F2_.point_reflection(center);
}

void Ellipse::reflect(const Line& axis) {
    F1_ = F1_.point_reflection(axis);
    F2_ = F2_.point_reflection(axis);
}

void Ellipse::scale(const Point& center, double coefficient) {
    F1_ = F1_.point_scale(center, coefficient);
    F2_ = F2_.point_scale(center, coefficient);
    a_ *= coefficient;
}

double Ellipse::get_small_axis() const {
    double c = F1_.distance_to(F2_) / 2;
    return sqrt(a_ * a_ - c * c);
}

bool Ellipse::containsPoint(const Point& point) {
    double curr_lenght = point.distance_to(F1_) + point.distance_to(F2_);
    const double EPS = 1e-9;
    return 2 * a_ - curr_lenght >= -EPS;
}

Polygon::Polygon(const std::vector<Point>& vertexes) {
    vertexes_ = vertexes;
}
size_t Polygon::verticesCount() const {
    return vertexes_.size();
}
std::vector<Point> Polygon::getVertices() const {
    return vertexes_;
}

bool Polygon::isConvex() {
    const double EPS = 1e-9;
    Point p1 = vertexes_[0], p2 = vertexes_[1], p3 = vertexes_[2];
    double curr = p1.x * p2.y + p2.x * p3.y + p3.x * p1.y - p1.y * p2.x -
                  p2.y * p3.x - p3.y * p1.x;
    curr /= fabs(curr);
    for (size_t i = 1; i < vertexes_.size(); ++i) {
        p1 = vertexes_[i], p2 = vertexes_[(i + 1) % vertexes_.size()];
        p3 = vertexes_[(i + 2) % vertexes_.size()];
        double temp = p1.x * p2.y + p2.x * p3.y + p3.x * p1.y - p1.y * p2.x -
                      p2.y * p3.x - p3.y * p1.x;
        if (curr * temp < -EPS) {
            return false;
        }
    }
    return true;
}

double Polygon::perimeter() const {
    double ans = 0;
    for (size_t i = 0; i < vertexes_.size(); ++i) {
        ans += vertexes_[i].distance_to(vertexes_[(i + 1) % vertexes_.size()]);
    }
    return ans;
}

double Polygon::area() const {
    return fabs(oriented_area());
}

double Polygon::oriented_area() const {
    double ans = 0;
    for (size_t i = 0; i < vertexes_.size(); ++i) {
        ans += vertexes_[i].x * vertexes_[(i + 1) % vertexes_.size()].y;
        ans -= vertexes_[i].y * vertexes_[(i + 1) % vertexes_.size()].x;
    }
    return ans / 2;
}

void Polygon::rotate(const Point& center, double angle) {
    for (size_t i = 0; i < vertexes_.size(); ++i) {
        vertexes_[i] = vertexes_[i].point_rotation(center, angle / 180 * M_PI);
    }
}
void Polygon::reflect(const Point& center) {
    for (size_t i = 0; i < vertexes_.size(); ++i) {
        vertexes_[i] = vertexes_[i].point_reflection(center);
    }
}

void Polygon::reflect(const Line& axis) {
    for (size_t i = 0; i < vertexes_.size(); ++i) {
        vertexes_[i] = vertexes_[i].point_reflection(axis);
    }
}

std::ostream& operator<<(std::ostream& os, const Polygon& pol) {
    std::vector<Point> v = pol.getVertices();
    for (size_t i = 0; i < v.size(); ++i) {
        os << v[i] << " ";
    }
    return os;
}

void Polygon::scale(const Point& center, double coefficient) {
    for (size_t i = 0; i < vertexes_.size(); ++i) {
        vertexes_[i] = vertexes_[i].point_scale(center, coefficient);
    }
}

bool Polygon::containsPoint(const Point& point) {
    double ans = 0;
    for (size_t i = 0; i < vertexes_.size(); ++i) {
        Point p1 = vertexes_[i], p2 = vertexes_[(i + 1) % vertexes_.size()],
              p3 = point;
        ans += fabs(p1.x * p2.y + p2.x * p3.y + p3.x * p1.y - p1.y * p2.x -
                    p2.y * p3.x - p3.y * p1.x) /
               2;
    }
    ans = fabs(ans);
    return equality(ans, area());
}

bool Rectangle::isConvex() {
    return true;
}

Rectangle::Rectangle(const Point& v1, const Point& v2, double k)
    : Polygon(v1, v2) {
    Line line_diagonal(v1, v2);
    double angle = (k < 1 ? atan(1 / k) : atan(k)) / M_PI * 180;
    Point v3 = line_intersections(
        rotation_line_about_point(line_diagonal, v1, angle / 180 * M_PI),
        rotation_line_about_point(line_diagonal, v2,
                                  (angle - 90) / 180 * M_PI));
    vertexes_.clear();
    vertexes_ = {v1, v3, v2};
    vertexes_.push_back(v3.point_reflection(center_of_segment(v1, v2)));
}

Point Rectangle::center() {
    return line_intersections(diagonals().first, diagonals().second);
}

std::pair<Line, Line> Rectangle::diagonals() {
    Line first_diag(vertexes_[0], vertexes_[2]);
    Line second_diag(vertexes_[1], vertexes_[3]);
    return {first_diag, second_diag};
}

bool Shape::operator!=(const Shape& another) const {
    return !(*this == another);
}

bool Shape::operator==(const Shape& another) const {
    const Polygon* this_ptr = dynamic_cast<const Polygon*>(this);
    const Polygon* another_ptr = dynamic_cast<const Polygon*>(&another);
    if (this_ptr == nullptr && another_ptr == nullptr) {
        const Ellipse* this_ptr = dynamic_cast<const Ellipse*>(this);
        const Ellipse* another_ptr = dynamic_cast<const Ellipse*>(&another);
        return ((this_ptr->focuses().first == another_ptr->focuses().first &&
                 this_ptr->focuses().second == another_ptr->focuses().second) ||
                (this_ptr->focuses().first == another_ptr->focuses().second &&
                 this_ptr->focuses().second == another_ptr->focuses().first)) &&
               equality(this_ptr->eccentricity(), another_ptr->eccentricity());
    }
    if (this_ptr != nullptr && another_ptr != nullptr) {
        if (this_ptr->verticesCount() != another_ptr->verticesCount()) {
            return false;
        }
        size_t size = this_ptr->verticesCount();
        std::vector<Point> this_vertices = this_ptr->getVertices();
        std::vector<Point> another_vertices = another_ptr->getVertices();
        for (size_t i = 0; i < size; ++i) {
            bool f = true;
            for (size_t j = 0; j < size; ++j) {
                if (this_vertices[j] != another_vertices[(j + i) % size]) {
                    f = false;
                    break;
                }
            }
            if (f) {
                return true;
            }
            f = true;
            for (size_t j = 0; j < size; ++j) {
                if (this_vertices[j] !=
                    another_vertices[(2 * size + 1 - j - i) % size]) {
                    f = false;
                    break;
                }
            }
            if (f) {
                return true;
            }
        }
        return false;
    }
    return false;
}

bool Shape::isCongruentTo(const Shape& another) {
    if (!equality(this->area(), another.area())) {
        return false;
    }
    return isSimilarTo(another);
}

bool Shape::isSimilarTo(const Shape& another) {
    const double EPS = 1e-9;
    Polygon* this_ptr = dynamic_cast<Polygon*>(this);
    const Polygon* another_ptr = dynamic_cast<const Polygon*>(&another);
    if (this_ptr == nullptr && another_ptr == nullptr) {
        Ellipse* this_ptr = dynamic_cast<Ellipse*>(this);
        const Ellipse* another_ptr = dynamic_cast<const Ellipse*>(&another);
        return fabs(this_ptr->eccentricity() - another_ptr->eccentricity()) <
               EPS;
    }
    if (this_ptr != nullptr && another_ptr != nullptr) {
        if (this_ptr->verticesCount() != another_ptr->verticesCount()) {
            return false;
        }
        double coeff = sqrt(this_ptr->area() / another_ptr->area());
        size_t size = this_ptr->verticesCount();
        std::vector<Point> this_vertices = this_ptr->getVertices();
        std::vector<Point> another_vertices = another_ptr->getVertices();
        for (size_t i = 0; i < size; ++i) {
            bool f = true;
            for (size_t j = 0; j < size; ++j) {
                double this_side =
                    this_vertices[j].distance_to(this_vertices[(j + 1) % size]);
                double another_side =
                    another_vertices[(i + j) % size].distance_to(
                        another_vertices[(i + j + 1) % size]);
                if (!(equality(this_side / another_side, coeff))) {
                    f = false;
                    break;
                }
            }
            if (f) {
                return true;
            }
            f = true;
            for (size_t j = 0; j < size; ++j) {
                double this_side =
                    this_vertices[j].distance_to(this_vertices[(j + 1) % size]);
                double another_side =
                    another_vertices[(2 * size + 1 - i - j) % size].distance_to(
                        another_vertices[(2 * size - i - j) % size]);
                if (!(equality(this_side / another_side, coeff))) {
                    f = false;
                    break;
                }
            }
            if (f) {
                return true;
            }
        }
        return false;
    }
    return false;
}

Circle Square::circumscribedCircle() {
    Point p_center = center();
    double radius = p_center.distance_to(vertexes_[0]);
    Circle ans(p_center, radius);
    return ans;
}

Circle Square::inscribedCircle() {
    Point p_center = center();
    double radius = vertexes_[1].distance_to(vertexes_[0]) / 2;
    Circle ans(p_center, radius);
    return ans;
}

Line Triangle::get_median(const Point& point) {
    Point center;
    if (point == vertexes_[0]) {
        center = center_of_segment(vertexes_[1], vertexes_[2]);
    } else if (point == vertexes_[1]) {
        center = center_of_segment(vertexes_[2], vertexes_[0]);
    } else {
        center = center_of_segment(vertexes_[0], vertexes_[1]);
    }
    return Line(center, point);
}

Circle Triangle::circumscribedCircle() {
    Circle ans(vertexes_[0], vertexes_[1], vertexes_[2]);
    return ans;
}

Circle Triangle::inscribedCircle() {
    Point center(line_intersections(get_bisector(vertexes_[0]),
                                    get_bisector(vertexes_[1])));
    double radius = center.distance_to(Line(vertexes_[0], vertexes_[1]));
    return Circle(center, 2 * radius);
}

Point Triangle::centroid() {
    Line first_median = get_median(vertexes_[0]);
    Line second_median = get_median(vertexes_[1]);
    return line_intersections(first_median, second_median);
}

Line Triangle::get_bisector(const Point& point) {
    double lenght1, lenght2, lymbda;
    Point base;
    if (point == vertexes_[0]) {
        lenght1 = vertexes_[0].distance_to(vertexes_[1]);
        lenght2 = vertexes_[0].distance_to(vertexes_[2]);
        lymbda = lenght1 / lenght2;
        base.x = (vertexes_[1].x + lymbda * vertexes_[2].x) / (1 + lymbda);
        base.y = (vertexes_[1].y + lymbda * vertexes_[2].y) / (1 + lymbda);
    } else if (point == vertexes_[1]) {
        lenght1 = vertexes_[1].distance_to(vertexes_[2]);
        lenght2 = vertexes_[1].distance_to(vertexes_[0]);
        lymbda = lenght1 / lenght2;
        base.x = (vertexes_[2].x + lymbda * vertexes_[0].x) / (1 + lymbda);
        base.y = (vertexes_[2].y + lymbda * vertexes_[0].y) / (1 + lymbda);
    } else {
        lenght1 = vertexes_[2].distance_to(vertexes_[0]);
        lenght2 = vertexes_[2].distance_to(vertexes_[1]);
        lymbda = lenght1 / lenght2;
        base.x = (vertexes_[0].x + lymbda * vertexes_[1].x) / (1 + lymbda);
        base.y = (vertexes_[0].y + lymbda * vertexes_[1].y) / (1 + lymbda);
    }
    return Line(point, base);
}

Point Triangle::orthocenter() {
    Line first_side(vertexes_[0], vertexes_[1]),
        second_side(vertexes_[1], vertexes_[2]);
    Line first_height = first_side.perpendicular_from_point(vertexes_[2]);
    Line second_hieght = second_side.perpendicular_from_point(vertexes_[0]);
    return line_intersections(first_height, second_hieght);
}

Line Triangle::EulerLine() {
    Point p_center = centroid();
    Point p_orthocenter = orthocenter();
    Line ans(p_center, p_orthocenter);
    return ans;
}

Circle Triangle::ninePointsCircle() {
    Point center_1 = center_of_segment(vertexes_[0], vertexes_[1]);
    Point center_2 = center_of_segment(vertexes_[1], vertexes_[2]);
    Point center_3 = center_of_segment(vertexes_[2], vertexes_[0]);
    return Circle(center_1, center_2, center_3);
}

bool Triangle::isConvex() {
    return true;
}
