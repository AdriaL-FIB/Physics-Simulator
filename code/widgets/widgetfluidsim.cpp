#include "widgetfluidsim.h"
#include "ui_widgetfluidsim.h"

WidgetFluidSim::WidgetFluidSim(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::WidgetFluidSim)
{
    ui->setupUi(this);

    connect(ui->btnUpdate, &QPushButton::clicked, this,
            [=] (void) { emit updatedParameters(); });
}

WidgetFluidSim::~WidgetFluidSim()
{
    delete ui;
}

double WidgetFluidSim::getGravity() const {
    return ui->gravity->value();
}

double WidgetFluidSim::getKBounce() const {
    return ui->kBounce->value();
}

double WidgetFluidSim::getKFriction() const {
    return ui->kFriction->value();
}

double WidgetFluidSim::getEmissionRate() const {
    return ui->emitRate->value();
}

double WidgetFluidSim::getMaxParticleLife() const {
    return ui->maxParticleLife->value();
}

double WidgetFluidSim::getParticleRadius() const {
    return ui->particleRadius->value();
}

double WidgetFluidSim::getParticleMass() const {
    return ui->particleMass->value();
}

bool WidgetFluidSim::emitParticles() const {
    return ui->emitParticles->isChecked();
}

double WidgetFluidSim::getRestDensity() const {
    return ui->restDensity->value();
}

double WidgetFluidSim::getViscosityConstant() const {
    return ui->mu->value();
}

int WidgetFluidSim::getIntegratorType() const {
    return ui->integrator->currentIndex();
}

void WidgetFluidSim::setSolverTypes(const std::vector<std::string>& solvers) {
    ui->integrator->clear();
    for (const std::string& s : solvers) {
        ui->integrator->addItem(QString::fromStdString(s));
    }
}