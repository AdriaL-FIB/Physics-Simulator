#ifndef WIDGETFOUNTAIN_H
#define WIDGETFOUNTAIN_H

#include <QWidget>

namespace Ui {
class WidgetFountain;
}

class WidgetFountain : public QWidget
{
    Q_OBJECT
public:
    explicit WidgetFountain(QWidget *parent = nullptr);
    ~WidgetFountain();

    double getGravity()    const;
    double getKBounce()    const;
    double getKFriction()    const;
    double getEmissionRate()    const;
    double getMaxParticleLife()    const;

    int    getIntegratorType() const;
    void setSolverTypes(const std::vector<std::string>& solvers);



signals:
    void updatedParameters();

private:
    Ui::WidgetFountain *ui;
};

#endif // WIDGETFOUNTAIN_H
