#include "widgetfountain.h"
#include "ui_widgetfountain.h"

WidgetFountain::WidgetFountain(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::WidgetFountain)
{
    ui->setupUi(this);

    connect(ui->btnUpdate, &QPushButton::clicked, this,
            [=] (void) { emit updatedParameters(); });
}

WidgetFountain::~WidgetFountain()
{
    delete ui;
}

double WidgetFountain::getGravity() const {
    return ui->gravity->value();
}

double WidgetFountain::getKBounce() const {
    return ui->kBounce->value();
}

double WidgetFountain::getKFriction() const {
    return ui->kFriction->value();
}

double WidgetFountain::getEmissionRate() const {
    return ui->emitRate->value();
}

double WidgetFountain::getMaxParticleLife() const {
    return ui->maxParticleLife->value();
}

int WidgetFountain::getIntegratorType() const {
    return ui->integrator->currentIndex();
}

void WidgetFountain::setSolverTypes(const std::vector<std::string>& solvers) {
    ui->integrator->clear();
    for (const std::string& s : solvers) {
        ui->integrator->addItem(QString::fromStdString(s));
    }
}