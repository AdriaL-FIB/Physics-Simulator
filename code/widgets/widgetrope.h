#ifndef WIDGETROPE_H
#define WIDGETROPE_H

#include <QWidget>

namespace Ui {
class WidgetRope;
}

class WidgetRope : public QWidget
{
    Q_OBJECT

public:
    explicit WidgetRope(QWidget *parent = nullptr);
    ~WidgetRope();

    double getGravity()       const;
    double getRopeLength()    const;
    double getStiffness()     const;
    double getDamping()       const;
    double getParticleScale() const;
    int getNumParticles()     const;
    bool checkCollisions()    const;
    bool anchorRopeEnd()      const;
    bool showParticles()      const;
    int getStartConfig()      const;
    int getSelectedIntegrator() const;

    void setSolverTypes(const std::vector<std::string>& solvers);

signals:
    void updatedParameters();

private:
    Ui::WidgetRope *ui;
};

#endif // WIDGETROPE_H
