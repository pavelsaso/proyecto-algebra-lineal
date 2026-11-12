#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <string>
#include <random>
using namespace Eigen;

struct Brazo2R {
    double l1, l2;

    // Cinematica directa: angulos -> posicion de la mano
    Vector2d fk(const Vector2d& q) const {
        return { l1*cos(q(0)) + l2*cos(q(0)+q(1)),
                 l1*sin(q(0)) + l2*sin(q(0)+q(1)) };
    }

    // Jacobiano: relaciona cambios de angulos con cambios de posicion
    Matrix2d jacobiano(const Vector2d& q) const {
        double s1 = sin(q(0)), c1 = cos(q(0));
        double s12 = sin(q(0)+q(1)), c12 = cos(q(0)+q(1));
        Matrix2d J;
        J << -l1*s1 - l2*s12, -l2*s12,
              l1*c1 + l2*c12,  l2*c12;
        return J;
    }

    // Metodo A: Gauss-Newton amortiguado con lambda fijo
    bool ikGN(const Vector2d& objetivo, Vector2d& q,
              double lambda = 0.01, int maxIter = 100, double tol = 1e-6) const {
        for (int i = 0; i < maxIter; ++i) {
            Vector2d e = objetivo - fk(q);
            if (e.norm() < tol) return true;
            Matrix2d J = jacobiano(q);
            Matrix2d A = J * J.transpose() + lambda*lambda*Matrix2d::Identity();
            q += J.transpose() * A.ldlt().solve(e);
        }
        return false;
    }

    // Metodo B: Levenberg-Marquardt (solo acepta pasos que reducen el error)
    bool ikLM(const Vector2d& objetivo, Vector2d& q,
              double lambda = 0.1, int maxIter = 200, double tol = 1e-6) const {
        double err = (objetivo - fk(q)).norm();
        for (int i = 0; i < maxIter; ++i) {
            if (err < tol) return true;
            Vector2d e = objetivo - fk(q);
            Matrix2d J = jacobiano(q);
            Matrix2d A = J * J.transpose() + lambda*lambda*Matrix2d::Identity();
            Vector2d qn = q + J.transpose() * A.ldlt().solve(e);
            double errN = (objetivo - fk(qn)).norm();
            if (errN < err) {
                q = qn; err = errN;
                lambda = std::max(lambda * 0.5, 1e-6);
            } else {
                lambda *= 2.0;
            }
        }
        return err < tol;
    }

    // IK analitica (ley de cosenos), rama codo abajo (q2 >= 0)
    bool ikAnalitica(const Vector2d& p, Vector2d& q) const {
        double c2 = (p.squaredNorm() - l1*l1 - l2*l2) / (2*l1*l2);
        if (std::abs(c2) > 1.0) return false;
        double q2 = std::acos(c2);
        double q1 = std::atan2(p(1), p(0))
                  - std::atan2(l2*sin(q2), l1 + l2*cos(q2));
        q << q1, q2;
        return true;
    }
};

// Prueba un metodo sobre 480 objetivos (15 radios x 32 angulos)
void barrido(const Brazo2R& brazo, bool usarLM, const std::string& nombre) {
    int total = 0, convergio = 0, codoAbajo = 0;
    double errMax = 0, errSuma = 0, tiempoSuma = 0;
    for (int i = 0; i < 15; ++i) {
        double r = 0.3 + 0.1*i;
        for (int j = 0; j < 32; ++j) {
            double th = 0.2*j;
            Vector2d p(r*cos(th), r*sin(th));
            Vector2d q(0.3, 0.5);   // punto inicial fijo
            auto t0 = std::chrono::high_resolution_clock::now();
            bool ok = usarLM ? brazo.ikLM(p, q) : brazo.ikGN(p, q);
            auto t1 = std::chrono::high_resolution_clock::now();
            tiempoSuma += std::chrono::duration<double>(t1 - t0).count();
            total++;
            if (ok) {
                convergio++;
                double err = (brazo.fk(q) - p).norm();
                errSuma += err;
                errMax = std::max(errMax, err);
                double q2 = std::atan2(std::sin(q(1)), std::cos(q(1)));
                if (q2 >= 0) codoAbajo++;
            }
        }
    }
    std::cout << "=== " << nombre << " ===\n";
    std::cout << "Objetivos probados : " << total << "\n";
    std::cout << "Convergieron       : " << convergio << "\n";
    std::cout << "En rama codo abajo : " << codoAbajo << "\n";
    std::cout << "Error medio (pos)  : " << errSuma / convergio << "\n";
    std::cout << "Error maximo (pos) : " << errMax << "\n";
    std::cout << "Tiempo medio por IK: " << tiempoSuma / total * 1e6 << " us\n\n";
}

// ===== Validaciones adicionales =====
static int convergencias(const Brazo2R& b, bool lm, double lam) {
    int c = 0;
    for (int i = 0; i < 15; ++i) for (int j = 0; j < 32; ++j) {
        double r = 0.3 + 0.1*i, th = 0.2*j;
        Vector2d p(r*cos(th), r*sin(th)), q(0.3, 0.5);
        if (lm ? b.ikLM(p, q, lam) : b.ikGN(p, q, lam)) c++;
    }
    return c;
}

void trayectoria2(const Brazo2R& b, double cx, double cy, double R, const char* arch) {
    const double PI = std::acos(-1.0);
    std::ofstream g(arch);
    g << "t,x_obj,y_obj,x,y,q1,q2,error,kappa\n";
    Vector2d q(0.3, 0.5); int fallos = 0; double eMax = 0, kMax = 0;
    for (int k = 0; k <= 200; ++k) {
        double t = 2*PI*k/200;
        Vector2d p(cx + R*cos(t), cy + R*sin(t));
        if (!b.ikLM(p, q)) fallos++;
        Vector2d pa = b.fk(q); double e = (pa - p).norm(); eMax = std::max(eMax, e);
        JacobiSVD<Matrix2d> s(b.jacobiano(q));
        double kap = s.singularValues()(0) / s.singularValues()(1); kMax = std::max(kMax, kap);
        g << t << "," << p(0) << "," << p(1) << "," << pa(0) << "," << pa(1) << ","
          << q(0) << "," << q(1) << "," << e << "," << kap << "\n";
    }
    std::cout << "Tray2 (" << cx << "," << cy << "," << R << "): fallos=" << fallos
              << " errmax=" << eMax << " kappa_max=" << kMax << "\n";
}

void validacionExtra(const Brazo2R& b) {
    const double PI = std::acos(-1.0);
    std::mt19937 rng(42);
    auto U = [&](double a, double c) { return a + (c - a) * (rng() / 4294967296.0); };
    std::cout << "\n##### VALIDACION ADICIONAL (semilla 42) #####\n";

    // A) Jacobiano analitico vs diferencias finitas centrales
    double maxJ = 0;
    for (int i = 0; i < 1000; ++i) {
        Vector2d q(U(-PI, PI), U(-PI, PI)); Matrix2d Jn; double h = 1e-6;
        for (int k = 0; k < 2; ++k) {
            Vector2d d = Vector2d::Zero(); d(k) = h;
            Jn.col(k) = (b.fk(q + d) - b.fk(q - d)) / (2*h);
        }
        maxJ = std::max(maxJ, (Jn - b.jacobiano(q)).cwiseAbs().maxCoeff());
    }
    std::cout << "A) max|J_analitico - J_dif.finitas| (1000 q) = " << maxJ << "\n";

    // B) Paso amortiguado: formula normal vs filtro SVD sigma/(sigma^2+lambda^2)
    double maxS = 0, lam = 0.1;
    for (int i = 0; i < 1000; ++i) {
        Vector2d q(U(-PI, PI), U(-PI, PI)), e(U(-1, 1), U(-1, 1));
        Matrix2d J = b.jacobiano(q);
        Vector2d d1 = J.transpose() * (J*J.transpose() + lam*lam*Matrix2d::Identity()).ldlt().solve(e);
        JacobiSVD<Matrix2d> s(J, ComputeFullU | ComputeFullV);
        Vector2d sv = s.singularValues();
        Vector2d f = (sv.array() / (sv.array().square() + lam*lam)).matrix();
        Vector2d d2 = s.matrixV() * (f.asDiagonal() * (s.matrixU().transpose() * e));
        maxS = std::max(maxS, (d1 - d2).norm());
    }
    std::cout << "B) max||paso_normal - paso_SVD|| (1000 casos) = " << maxS << "\n";

    // C) Robustez: 1000 objetivos aleatorios y puntos iniciales aleatorios
    int okGN = 0, okLM = 0;
    for (int i = 0; i < 1000; ++i) {
        double r = U(0.25, 1.75), th = U(-PI, PI);
        Vector2d p(r*cos(th), r*sin(th)), q0(U(-PI, PI), U(-PI, PI)), q = q0;
        if (b.ikGN(p, q)) okGN++;
        q = q0; if (b.ikLM(p, q)) okLM++;
    }
    std::cout << "C) Inicio aleatorio, 1000 casos: GN=" << okGN << " LM=" << okLM << "\n";

    // D) Sensibilidad al parametro lambda (malla de 480 objetivos, inicio fijo)
    std::cout << "D) lambda  GN  LM\n";
    for (double l : {0.001, 0.01, 0.1, 0.5, 1.0})
        std::cout << "D) " << l << " " << convergencias(b, false, l) << " " << convergencias(b, true, l) << "\n";

    // E) Amplificacion de error: ||dq||/||dp|| vs 1/sigma_min(J)
    std::cout << "E) q2  sigma_min  1/sigma_min  max_ratio_medido  kappa\n";
    for (double q2 : {2.35, 1.5, 0.8, 0.3, 0.1}) {
        Vector2d q(0.5, q2), p = b.fk(q), q0;
        b.ikAnalitica(p, q0);
        JacobiSVD<Matrix2d> s(b.jacobiano(q0));
        double eps = 1e-7, rmax = 0;
        for (int k = 0; k < 360; ++k) {
            double ph = 2*PI*k/360; Vector2d d(eps*cos(ph), eps*sin(ph)), q1;
            b.ikAnalitica(p + d, q1);
            rmax = std::max(rmax, (q1 - q0).norm() / eps);
        }
        std::cout << "E) " << q2 << " " << s.singularValues()(1) << " " << 1/s.singularValues()(1)
                  << " " << rmax << " " << s.singularValues()(0)/s.singularValues()(1) << "\n";
    }

    // F) Mapa de convergencia (para figura) y segunda trayectoria
    std::ofstream m("mapa.csv"); m << "x,y,r,gn,lm\n";
    for (int i = 0; i < 15; ++i) for (int j = 0; j < 32; ++j) {
        double r = 0.3 + 0.1*i, th = 0.2*j;
        Vector2d p(r*cos(th), r*sin(th)), q(0.3, 0.5), q2(0.3, 0.5);
        int g = b.ikGN(p, q), l = b.ikLM(p, q2);
        m << p(0) << "," << p(1) << "," << r << "," << g << "," << l << "\n";
    }
    trayectoria2(b, 1.2, 0.0, 0.55, "trayectoria2.csv");
}

int main() {
    Brazo2R brazo{1.0, 0.8};
    const double PI = std::acos(-1.0);

    // ---- 1) Un objetivo: numerica vs analitica ----
    Vector2d obj(1.2, 0.9), qn(0.3, 0.5), qa;
    brazo.ikLM(obj, qn);
    brazo.ikAnalitica(obj, qa);
    std::cout << "q numerica  = " << qn.transpose() << "\n";
    std::cout << "q analitica = " << qa.transpose() << "\n";
    std::cout << "||q_num - q_ana|| = " << (qn - qa).norm() << "\n\n";

    // ---- 2) Comparacion de metodos sobre 480 objetivos ----
    barrido(brazo, false, "Gauss-Newton amortiguado (lambda fijo)");
    barrido(brazo, true,  "Levenberg-Marquardt adaptativo");

    // ---- 3) Barrido de kappa(J) vs q2 -> CSV ----
    std::ofstream f("kappa.csv");
    f << "q2,det,sigma_max,sigma_min,kappa\n";
    for (double q2 = -PI; q2 <= PI; q2 += 0.01) {
        Vector2d q(0.0, q2);
        Matrix2d J = brazo.jacobiano(q);
        JacobiSVD<Matrix2d> svd(J);
        double smax = svd.singularValues()(0), smin = svd.singularValues()(1);
        double kappa = (smin < 1e-12) ? 1e12 : smax / smin;
        f << q2 << "," << J.determinant() << "," << smax << ","
          << smin << "," << kappa << "\n";
    }
    std::cout << "kappa.csv generado\n";

    // ---- 4) Trayectoria circular (cada IK parte de la solucion anterior) ----
    std::ofstream g("trayectoria.csv");
    g << "t,x_obj,y_obj,x,y,q1,q2,error,kappa\n";
    const int N = 200;
    const double cx = 1.0, cy = 0.5, R = 0.4;
    Vector2d q(0.3, 0.5);
    int fallos = 0;
    double errMaxT = 0;
    for (int k = 0; k <= N; ++k) {
        double t = 2*PI*k/N;
        Vector2d p(cx + R*cos(t), cy + R*sin(t));
        bool ok = brazo.ikLM(p, q);
        if (!ok) fallos++;
        Vector2d pa = brazo.fk(q);
        double err = (pa - p).norm();
        errMaxT = std::max(errMaxT, err);
        JacobiSVD<Matrix2d> s(brazo.jacobiano(q));
        double kap = s.singularValues()(0) / s.singularValues()(1);
        g << t << "," << p(0) << "," << p(1) << "," << pa(0) << "," << pa(1)
          << "," << q(0) << "," << q(1) << "," << err << "," << kap << "\n";
    }
    std::cout << "Trayectoria: " << N+1 << " puntos, fallos = " << fallos
              << ", error max = " << errMaxT << "\ntrayectoria.csv generado\n";

    validacionExtra(brazo);
}