#ifndef WIDGETCLOTH_H
#define WIDGETCLOTH_H

#include <QWidget>
#include "defines.h"

namespace Ui {
class WidgetCloth;
}

class WidgetCloth : public QWidget
{
    Q_OBJECT

public:
    explicit WidgetCloth(QWidget *parent = nullptr);
    ~WidgetCloth();

    double getGravity()        const;
    double getStiffness()      const;
    double getDamping()        const;

    Vec2 getDimensions()       const;
    Vec2i getNumParticles()    const;
    double getParticleRadius() const;
    int getRelaxationSteps()   const;

    bool showParticles()       const;

    int getIntegratorType() const;
    void setSolverTypes(const std::vector<std::string>& solvers);

signals:
    void updatedParameters();
    void freeAnchors();

private:
    Ui::WidgetCloth *ui;
};

#endif // WIDGETCLOTH_H
