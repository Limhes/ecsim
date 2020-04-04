#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTableView>
#include <QStyledItemDelegate>
#include <QItemDelegate>

#include "simulation.h"
#include "qcustomplot.h"
#include "model.h"
#include "delegates.h"

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    void drawMetricsLine(bool, double, QString, QColor);
    void showError(QString, QString);

private slots:
    void removeSpecies();
    void removeRedox();
    void removeReaction();
    void cmbElectrodeTypeChanged(int);
    void edtVerticesChanged();
    void runSimulation();
    void removeAllPlots();
    void removeSelectedPlots();
    void addTheoreticalValuesToPlot();

private:
    void attachModelAndDelegateToView();
    void addPlot(vector<double>&, vector<double>&);
    void saveSystemToFile();
    void loadSystemFromFile();

    Simulation *sim;
    Model *model;

    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
