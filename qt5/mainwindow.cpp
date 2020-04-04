#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QPushButton>
#include <QMessageBox>
#include <QDebug>
#include <QVector>
#include <cmath>

MainWindow::MainWindow(QWidget *parent)
  : QMainWindow(parent), ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    sim = nullptr;
    model = new Model();
    model->loadDefaultValues();

    attachModelAndDelegateToView();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::attachModelAndDelegateToView()
{
    // delegate for species:
    QStyledItemDelegate *delegateSpeciesItem = new SpeciesItemDelegate(this, model->sys);
    // delegate for enabled checkbox:
    QItemDelegate *delegateCheckBoxItem = new CheckBoxItemDelegate(this);
    // delegate for rate constants (high precision unsigned double):
    QStyledItemDelegate *delegateRateConstantItem = new QStyledItemDelegate(this);
    DoubleFactory *factoryRateConstantItem = new DoubleFactory();
    factoryRateConstantItem->setScientificFormat(true, 1.0e-6); // number is the smallest rate constant allowed (apart from 0)
    factoryRateConstantItem->setSpinBoxLimits(0.0, 1.0e10, 3);
    delegateRateConstantItem->setItemEditorFactory(factoryRateConstantItem);
    // delegate for potentials (low precision signed double):
    QStyledItemDelegate *delegatePotentialItem = new QStyledItemDelegate(this);
    DoubleFactory *factoryPotentialItem = new DoubleFactory();
    factoryPotentialItem->setSpinBoxLimits(-10.0, 10.0, 3);
    delegatePotentialItem->setItemEditorFactory(factoryPotentialItem);
    // delegate for number of electrons (natural number, not too big):
    QStyledItemDelegate *delegateNumElectronsItem = new QStyledItemDelegate(this);
    IntFactory *factoryNumElectronsItem = new IntFactory();
    factoryNumElectronsItem->setSpinBoxLimits(1, 100);
    delegateNumElectronsItem->setItemEditorFactory(factoryNumElectronsItem);
    // delegate for alpha (double between 0 and 1):
    QStyledItemDelegate *delegateAlphaItem = new QStyledItemDelegate(this);
    DoubleFactory *factoryAlphaItem = new DoubleFactory();
    factoryAlphaItem->setSpinBoxLimits(0.0, 1.0, 3);
    delegateAlphaItem->setItemEditorFactory(factoryAlphaItem);
    // delegate for concentration [mol/m3] (double between 0 and 100000):
    QStyledItemDelegate *delegateConcentrationItem = new QStyledItemDelegate(this);
    DoubleFactory *factoryConcentrationItem = new DoubleFactory();
    factoryConcentrationItem->setSpinBoxLimits(0.0, 1.0e5, 3);
    delegateConcentrationItem->setItemEditorFactory(factoryConcentrationItem);
    // delegate for diffusion coefficient [m2/s] (double between 0 and 100000):
    QStyledItemDelegate *delegateDiffCoeffItem = new QStyledItemDelegate(this);
    DoubleFactory *factoryDiffCoeffItem = new DoubleFactory();
    factoryDiffCoeffItem->setScientificFormat(true, 1.0e-12); // number is the smallest diff coeff allowed
    factoryDiffCoeffItem->setSpinBoxLimits(1.0e-12, 1.0e-6, 3);
    delegateDiffCoeffItem->setItemEditorFactory(factoryDiffCoeffItem);

    // redox table/buttons:
    ui->tblRedox->setModel(model->tblRedox);
    ui->tblRedox->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    ui->tblRedox->setItemDelegateForColumn(0, delegateCheckBoxItem);
    ui->tblRedox->setItemDelegateForColumn(1, delegateSpeciesItem);
    ui->tblRedox->setItemDelegateForColumn(2, delegateSpeciesItem);
    ui->tblRedox->setItemDelegateForColumn(3, delegatePotentialItem);
    ui->tblRedox->setItemDelegateForColumn(4, delegateNumElectronsItem);
    ui->tblRedox->setItemDelegateForColumn(5, delegateRateConstantItem);
    ui->tblRedox->setItemDelegateForColumn(6, delegateAlphaItem);
    connect(ui->btnAddRedox, SIGNAL(clicked()), model->tblRedox, SLOT(addData()));
    connect(ui->btnRemoveRedox, SIGNAL(clicked()), this, SLOT(removeRedox()));

    // reaction table/buttons:
    ui->tblReaction->setModel(model->tblReaction);
    ui->tblReaction->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    ui->tblReaction->setItemDelegateForColumn(0, delegateCheckBoxItem);
    ui->tblReaction->setItemDelegateForColumn(1, delegateSpeciesItem);
    ui->tblReaction->setItemDelegateForColumn(2, delegateSpeciesItem);
    ui->tblReaction->setItemDelegateForColumn(3, delegateSpeciesItem);
    ui->tblReaction->setItemDelegateForColumn(4, delegateSpeciesItem);
    ui->tblReaction->setItemDelegateForColumn(5, delegateRateConstantItem);
    ui->tblReaction->setItemDelegateForColumn(6, delegateRateConstantItem);
    connect(ui->btnAddReaction, SIGNAL(clicked()), model->tblReaction, SLOT(addData()));
    connect(ui->btnRemoveReaction, SIGNAL(clicked()), this, SLOT(removeReaction()));

    // species table/buttons:
    ui->tblSpecies->setModel(model->tblSpecies);
    ui->tblSpecies->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    // for column 0, create delegate that does not allow to pick a species name that already exists...
    ui->tblSpecies->setItemDelegateForColumn(1, delegateConcentrationItem);
    ui->tblSpecies->setItemDelegateForColumn(2, delegateDiffCoeffItem);
    connect(ui->btnAddSpecies, SIGNAL(clicked()), model->tblSpecies, SLOT(addData()));
    connect(ui->btnRemoveSpecies, SIGNAL(clicked()), this, SLOT(removeSpecies()));

    // update redox and reaction table views when a species name changes:
    connect(model->tblSpecies, SIGNAL(speciesNameChanged()), ui->tblRedox->viewport(), SLOT(update()));
    connect(model->tblSpecies, SIGNAL(speciesNameChanged()), ui->tblReaction->viewport(), SLOT(update()));

    // add electrode types:
    for (string s: Electrode::electrodeTypes)
        ui->cmbElectrodeType->addItem(QString::fromStdString(s));
    connect(ui->cmbElectrodeType, SIGNAL(currentIndexChanged(int)), this, SLOT(cmbElectrodeTypeChanged(int)));
    connect(ui->dblElectrodeGeometry1, SIGNAL(valueChanged(double)), model->el, SLOT(setGeom1(double)));
    connect(ui->dblElectrodeGeometry2, SIGNAL(valueChanged(double)), model->el, SLOT(setGeom2(double)));
    ui->dblElectrodeGeometry1->setValue(model->el->m_geom1);
    ui->dblElectrodeGeometry2->setValue(model->el->m_geom2);
    emit ui->cmbElectrodeType->currentIndexChanged(static_cast<int>(model->el->m_type));

    // add environment to GUI:
    connect(ui->dblTemperature, SIGNAL(valueChanged(double)), model->env, SLOT(setTemperature(double)));
    ui->dblTemperature->setValue(model->env->temperature);

    // add experiment to GUI:
    connect(ui->edtVertices, SIGNAL(editingFinished()), this, SLOT(edtVerticesChanged()));
    connect(ui->dblScanRate, SIGNAL(valueChanged(double)), model->exper, SLOT(setScanRate(double)));
    connect(ui->dblInitialPotential, SIGNAL(valueChanged(double)), model->exper, SLOT(setInitialPotential(double)));
    connect(ui->dblFinalPotential, SIGNAL(valueChanged(double)), model->exper, SLOT(setFinalPotential(double)));
    connect(ui->intNumCycles, SIGNAL(valueChanged(int)), model->exper, SLOT(setNumCycles(int)));
    ui->dblScanRate->setValue(model->exper->scanRate);
    ui->dblInitialPotential->setValue(model->exper->initialPotential);
    ui->dblFinalPotential->setValue(model->exper->finalPotential);
    ui->intNumCycles->setValue(model->exper->numCycles);
    ui->edtVertices->setText("-1.0");
    edtVerticesChanged();

    // button to run simulation and remove curves from plot:
    connect(ui->btnSimulate, SIGNAL(clicked()), this, SLOT(runSimulation()));
    connect(ui->btnRemoveAllPlots, SIGNAL(clicked()), this, SLOT(removeAllPlots()));
    connect(ui->btnRemoveSelectedPlots, SIGNAL(clicked()), this, SLOT(removeSelectedPlots()));
    connect(ui->btnAddPeakCurrentPotential, SIGNAL(clicked()), this, SLOT(addTheoreticalValuesToPlot()));

    // QCustomPlot stuff:
    ui->plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->plot->axisRect()->setupFullAxesBox();
    ui->plot->legend->setVisible(true);
    ui->plot->legend->setFont(QFont("Helvetica", 9));
    ui->plot->xAxis->setLabel("Potential [V]");
    ui->plot->yAxis->setLabel("Current [A] / sqrt(scan rate [V/s])");
}

void MainWindow::removeSelectedPlots()
{
    for (int i = 0; i < ui->plot->plottableCount(); ++i)
    {
        QCPAbstractPlottable *p = ui->plot->plottable(i);
        if (p->selected())
        {
            ui->plot->removePlottable(p);
        }
    }
    ui->plot->replot();
}

void MainWindow::removeAllPlots()
{
    ui->plot->clearPlottables(); // curves (CVs)
    ui->plot->clearItems(); // also clear all lines

    ui->plot->replot();
}

void MainWindow::drawMetricsLine(bool directionIsHorizontal, double value, QString str, QColor c)
{
    double imin, imax, Emin, Emax, itext, Etext;
    if (directionIsHorizontal) // plot a horizontal line
    {
        imin = value;
        imax = value;
        itext = value;
        Emin = ui->plot->xAxis->range().lower;
        Emax = ui->plot->xAxis->range().upper;
        Etext = Emin+(Emax-Emin)*0.05;
    }
    else // plot a vertical line
    {
        imin = ui->plot->yAxis->range().lower;
        imax = ui->plot->yAxis->range().upper;
        itext = imin+(imax-imin)*0.05;
        Emin = value;
        Emax = value;
        Etext = value;
    }

    QCPItemLine *p = new QCPItemLine(ui->plot);
    QPen pen(c);
    pen.setStyle(Qt::PenStyle::DashLine);
    p->setPen(pen);
    p->setSelectable(false);
    p->start->setType(QCPItemPosition::ptPlotCoords);
    p->end->setType(QCPItemPosition::ptPlotCoords);
    p->start->setCoords(Emin, imin);
    p->end->setCoords(Emax, imax);
    p->setTail(QCPLineEnding::esNone);
    p->setHead(QCPLineEnding::esNone);

    QCPItemText *t = new QCPItemText(ui->plot);
    t->setPen(Qt::NoPen);
    t->setColor(c);
    t->setBrush(QBrush(QColor("white")));
    t->setSelectable(false);
    t->position->setType(QCPItemPosition::ptPlotCoords);
    t->position->setCoords(Etext, itext);
    t->setText(str);
}

void MainWindow::addTheoreticalValuesToPlot()
{
    // add calculated reversible I-V locations as dots:
    double iPeak, EPeak;

    // metrics on first redox reaction:
    const double f = CONST_F / (CONST_R * model->env->temperature);
    const double iNorm = CONST_F*model->el->electrodeArea*sqrt(f);
    const double nu = model->exper->scanRate;
    const double scanDirection = (model->exper->vertexPotentials[0] - model->exper->initialPotential) / fabs(model->exper->vertexPotentials[0] - model->exper->initialPotential);
    Redox *redox = model->sys->vecAllRedox[0];
    const double conc = redox->specOxidized->initialConcentration;
    const double DO = redox->specOxidized->diffusionConstant;
    const double DR = redox->specReduced->diffusionConstant;
    const double k0 = redox->rateConstantHetero;
    const double alpha = redox->alpha;
    const double n = redox->numberElectrons;
    const double E0 = redox->standardPotential;
    const double Ehalf = E0 + 1/(2*n*f)*log(DR/DO);
    const double Lambda = redox->rateConstantHetero / sqrt( pow(DO, 1-redox->alpha)*pow(DR, redox->alpha)*f*nu );
    // dimensionless i/E:
    const double currentDO = 0.4463; // "Randles-Sevcik" current
    const double potentialDO = 1.109;
    const double currentKP = 0.4958; // Nicholson 1964
    const double potentialKP = 0.780;

    // remove all previous metrics info:
    ui->plot->clearItems();

    // Er:
    EPeak = Ehalf + scanDirection*potentialDO/(f*n);
    iPeak = scanDirection*currentDO*iNorm*conc*pow(n, 1.5);
    if (scanDirection < 0) iPeak *= sqrt(DO); else iPeak *= sqrt(DR);
    drawMetricsLine(true, iPeak, "i rev", QColor("black"));
    drawMetricsLine(false, EPeak, "E rev", QColor("black"));

    // Eirr:
    if (Lambda < pow(10.0, -2*(1+alpha)))
    {
        double D, beta;
        if (scanDirection < 0) { D = DO; beta = alpha; } else { D = DR; beta = 1-alpha; }
        iPeak = scanDirection*currentKP*iNorm*conc*sqrt(alpha)*sqrt(D);
        EPeak = E0 + scanDirection/f * ( potentialKP/alpha - log(k0 * sqrt(1/(beta*f*nu*D)))/beta ); // Saveant, Elements of ... (page 53)
        drawMetricsLine(true, iPeak, "i irrev", QColor("blue"));
        drawMetricsLine(false, EPeak, "E irrev", QColor("blue"));
    }

    if (model->sys->vecAllReactions.size() > 0)
    {
        Reaction *rxn = model->sys->vecAllReactions[0];
        const double lambda = (rxn->rateConstantForward + rxn->rateConstantBackward)/(f*nu);

        if (lambda > 5)
        {
            EPeak = Ehalf + scanDirection/f * ( potentialKP - 0.5*log(rxn->rateConstantForward/(f*nu)) ); // ErCi in KP: Saveant, Elements of ... (page 83)
            iPeak = scanDirection*currentKP*iNorm*conc;
            if (scanDirection < 0) iPeak *= sqrt(DO); else iPeak *= sqrt(DR);
            drawMetricsLine(true, iPeak, "i ErCi in KP", QColor("red"));
            drawMetricsLine(false, EPeak, "E ErCi in KP", QColor("red"));

            if (rxn->rateConstantBackward > MIN_RATE)
            {
                double K = rxn->rateConstantForward / rxn->rateConstantBackward;
                EPeak = Ehalf + scanDirection/f*(potentialDO - log(1+K)); // ErCr in DE: Saveant, Elements of ... (page 85)
                iPeak = scanDirection*currentDO*iNorm*conc;
                if (scanDirection < 0) iPeak *= sqrt(DO); else iPeak *= sqrt(DR);
                drawMetricsLine(true, iPeak, "i ErCr in DE", QColor("green"));
                drawMetricsLine(false, EPeak, "E ErCr in DE", QColor("green"));
            }
        }
    }

    // refresh plot:
    ui->plot->replot();
}

void MainWindow::addPlot(vector<double> &_current, vector<double> &_potential)
{
    // continuous color selection:
    QColor colours[8] = {QColor("magenta"), QColor("red"),
                          QColor("darkRed"), QColor("darkCyan"), QColor("darkMagenta"),
                          QColor("green"), QColor("darkGreen"), QColor("blue")};
    static size_t plotidx = 0;
    QColor color = colours[plotidx % 8];
    plotidx++;

    static size_t plottableidx = 0;
    plottableidx++;

    // display curve on plot:
    QCPCurve *cv = new QCPCurve(ui->plot->xAxis, ui->plot->yAxis);
    QVector<QCPCurveData> dataCV(static_cast<int>(_current.size()));
    for (size_t i = 0; i < _current.size(); i++)
    {
        dataCV[static_cast<int>(i)] = QCPCurveData(i, _potential[i], _current[i]);
    }
    cv->data()->set(dataCV, true);
    cv->setPen(QPen(color));
    cv->setName("Simulation ..."); // add plotidx into string

    // let the ranges scale themselves so graph 0 fits perfectly in the visible area:
    ui->plot->rescaleAxes();
    // zoom out a bit:
    ui->plot->yAxis->scaleRange(1.2, ui->plot->yAxis->range().center());
    ui->plot->xAxis->scaleRange(1.2, ui->plot->xAxis->range().center());
    // refresh plot:
    ui->plot->replot();
}

void MainWindow::showError(QString title, QString message)
{
    QMessageBox::critical(this, title, message, QMessageBox::Ok);
}

void MainWindow::runSimulation()
{
    vector<double> current, potential;
    stringstream simulationOutput(ios::out);

    // run checks:
    QString titleCannotRunSim = "Cannot run simulation";
    if (model->el->electrodeArea < MIN_AREA)
    {
        showError(titleCannotRunSim, "The electrode area is too small.");
        return;
    }
    size_t redoxEnabledCount = 0;
    for (auto redox: model->sys->vecAllRedox)
    {
        if (redox->enabled)
        {
            redoxEnabledCount++;
            if (redox->specReduced == nullptr || redox->specOxidized == nullptr)
            {
                showError(titleCannotRunSim, "Each (enabled!) redox step needs both Ox and Red species.");
                return;
            }
        }
    }
    if (redoxEnabledCount == 0)
    {
        showError(titleCannotRunSim, "The simulation needs at least one (enabled!) redox step.");
        return;
    }
    for (auto rxn: model->sys->vecAllReactions)
    {
        if (rxn->enabled && ((rxn->specLHS1 == nullptr && rxn->specLHS2 == nullptr) || (rxn->specRHS1 == nullptr && rxn->specRHS2 == nullptr)))
        {
            showError(titleCannotRunSim, "Each (enabled!) homogeneous reaction needs at least one species at either side of the arrow.");
            return;
        }
    }

    delete sim;
    // setup new simulation using the objects stored in the model object:
    sim = new Simulation(simulationOutput, model->env, model->el, model->exper, model->sys);
    sim->setGridSizing(ui->dblSimGamma->value(), ui->dblSimMinF->value(), ui->dblSimMaxF->value(),
                       ui->dblSimMinLogRate->value(), ui->dblSimMaxLogRate->value());
    sim->setPotentialSizing(ui->dblSimDeltaTheta->value());
    sim->setDifferentialOrders(static_cast<size_t>(ui->intSimNumDeriv->value()),
                               static_cast<size_t>(ui->intSimNumCurrent->value()));
    // run the simulation:
    sim->run(current, potential);

    // display simulation output:
    ui->txtSimulationOutput->setHtml("<html><head/><body><p>" + QString::fromStdString(simulationOutput.str()) + "</p></body></html>");

    // scale current by sqrt(scan rate):
    for (size_t x = 0; x < current.size(); x++) current[x] /= sqrt(model->exper->scanRate);
    // display plot:
    addPlot(current, potential);
}

void MainWindow::edtVerticesChanged()
{
    bool inputIsValid;
    double val;
    QString formattedVertices("");

    model->exper->vertexPotentials.clear();

    QStringList strlist = ui->edtVertices->text().split(" ", QString::SkipEmptyParts);
    for (QString part: strlist)
    {
        val = part.toDouble(&inputIsValid);
        if (inputIsValid)
        {
            model->exper->vertexPotentials.emplace_back(val);
            formattedVertices += QString{"%1 "}.arg(val, 0, 'f', 3);
        }
    }
    formattedVertices = formattedVertices.trimmed();

    ui->edtVertices->setText(formattedVertices);
}

void MainWindow::cmbElectrodeTypeChanged(int e)
{
    size_t elidx = static_cast<size_t>(e);

    model->el->setType(elidx); // set type and update()

    ui->lblElectrodeGeometry1->setText(QString::fromStdString(Electrode::electrodeGeom1[elidx]));
    ui->lblElectrodeGeometry2->setText(QString::fromStdString(Electrode::electrodeGeom2[elidx]));
    ui->dblElectrodeGeometry2->setDisabled((Electrode::electrodeGeom2[elidx] == NO_GEOM));
}

void MainWindow::removeSpecies()
{
    bool speciesInUse = false;
    QItemSelectionModel *select = ui->tblSpecies->selectionModel();

    if (select->hasSelection()) //check if has selection
    {
        vector<bool> selectedRows(model->sys->vecAllSpecies.size(), false);

        for (auto selectedIndexes: select->selectedIndexes())
        {
            selectedRows[static_cast<size_t>(selectedIndexes.row())] = true;
        }

        for (size_t s = 0; s < model->sys->vecAllSpecies.size(); s++)
        {
            auto spec = model->sys->vecAllSpecies[s];
            if (selectedRows[s]) // species is selected to be removed
            {
                // check if species is in use
                for (auto redox: model->sys->vecAllRedox)
                {
                    if (redox->specReduced == spec || redox->specOxidized == spec)
                        speciesInUse = true; // species is in use in a redox step
                }
                for (auto rxn: model->sys->vecAllReactions)
                {
                    if (rxn->specLHS1 == spec || rxn->specLHS2 == spec || rxn->specRHS1 == spec || rxn->specRHS2 == spec)
                        speciesInUse = true; // species is in use in a reaction
                }
            }
        }

        if (speciesInUse)
        {
            showError("Species in use", "One or more of the selected species are in use in a redox step or reaction. Remove species from redox step or reaction, even if they are disabled.");
            // show dialog that species is in use
        }
        else
        {
            for (int s = static_cast<int>(model->sys->vecAllSpecies.size())-1; s >= 0; s--) // delete from the bottom upwards!
            {
                if (selectedRows[static_cast<size_t>(s)]) // species is selected to be removed
                {
                    qDebug() << "deleting species at row " << s;
                    model->tblSpecies->removeRow(s);
                }
            }
        }
    }
}

void MainWindow::removeRedox()
{
    QItemSelectionModel *select = ui->tblRedox->selectionModel();

    if (select->hasSelection()) //check if has selection
    {
        vector<bool> selectedRows(model->sys->vecAllRedox.size(), false);

        for (auto selectedIndexes: select->selectedIndexes())
        {
            selectedRows[static_cast<size_t>(selectedIndexes.row())] = true;
        }

        for (int s = static_cast<int>(model->sys->vecAllRedox.size())-1; s >= 0; s--) // delete from the bottom upwards!
        {
            if (selectedRows[static_cast<size_t>(s)]) // redox is selected to be removed
            {
                qDebug() << "deleting redox at row " << s;
                model->tblRedox->removeRow(s);
            }
        }
    }
}

void MainWindow::removeReaction()
{
    QItemSelectionModel *select = ui->tblReaction->selectionModel();

    if (select->hasSelection()) //check if has selection
    {
        vector<bool> selectedRows(model->sys->vecAllReactions.size(), false);

        for (auto selectedIndexes: select->selectedIndexes())
        {
            selectedRows[static_cast<size_t>(selectedIndexes.row())] = true;
        }

        for (int s = static_cast<int>(model->sys->vecAllReactions.size())-1; s >= 0; s--) // delete from the bottom upwards!
        {
            if (selectedRows[static_cast<size_t>(s)]) // reaction is selected to be removed
            {
                qDebug() << "deleting reaction at row " << s;
                model->tblReaction->removeRow(s);
            }
        }
    }
}

void MainWindow::saveSystemToFile()
{

}
void MainWindow::loadSystemFromFile()
{

}
