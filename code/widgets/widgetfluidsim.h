#ifndef WIDGETFLUIDSIM_H
#define WIDGETFLUIDSIM_H

#include <QWidget>

namespace Ui {
class WidgetFluidSim;
}

class WidgetFluidSim : public QWidget
{
    Q_OBJECT
public:
    explicit WidgetFluidSim(QWidget *parent = nullptr);
    ~WidgetFluidSim();

    double getGravity()    const;
    double getKBounce()    const;
    double getKFriction()    const;
    double getEmissionRate()    const;
    double getMaxParticleLife()    const;
    double getParticleRadius()    const;
    double getParticleMass()    const;
    bool emitParticles()    const;
    double getRestDensity()    const;
    double getViscosityConstant()    const;
    bool getEnableWind()    const;

    int    getIntegratorType() const;
    void   setIntegratorType(int i);
    void setSolverTypes(const std::vector<std::string>& solvers);



signals:
    void updatedParameters();

private:
    Ui::WidgetFluidSim *ui;
};

#endif // WIDGETFLUIDSIM_H
