#ifndef MODEL_H
#define MODEL_H

#include <QStandardItemModel>
#include "simulation.h"

class TableModelRedox : public QAbstractTableModel
{
    Q_OBJECT

public:
    explicit TableModelRedox(QObject *parent = nullptr, System* sys = nullptr)
        : QAbstractTableModel(parent), m_data(sys) { }

    // Basic QAbstractTableModel API
    int rowCount(const QModelIndex&) const override;
    int columnCount(const QModelIndex&) const override;
    QVariant data(const QModelIndex&, int) const override;
    bool setData(const QModelIndex&, const QVariant&, int) override;
    Qt::ItemFlags flags(const QModelIndex&) const override;
    QVariant headerData(int, Qt::Orientation, int) const override;

    bool insertRows(int row, int count, const QModelIndex &parent) override;
    bool removeRows(int row, int count, const QModelIndex &parent) override;

public slots:
    void addData();

private:
    System* m_data;
};

class TableModelReaction : public QAbstractTableModel
{
    Q_OBJECT

public:
    explicit TableModelReaction(QObject *parent = nullptr, System* sys = nullptr)
        : QAbstractTableModel(parent), m_data(sys) { }

    // Basic QAbstractTableModel API
    int rowCount(const QModelIndex&) const override;
    int columnCount(const QModelIndex&) const override;
    QVariant data(const QModelIndex&, int) const override;
    bool setData(const QModelIndex&, const QVariant&, int) override;
    Qt::ItemFlags flags(const QModelIndex&) const override;
    QVariant headerData(int, Qt::Orientation, int) const override;

    bool insertRows(int row, int count, const QModelIndex &parent) override;
    bool removeRows(int row, int count, const QModelIndex &parent) override;

public slots:
    void addData();
private:
    System* m_data;
};

class TableModelSpecies : public QAbstractTableModel
{
    Q_OBJECT

public:
    explicit TableModelSpecies(QObject *parent = nullptr, System* sys = nullptr)
        : QAbstractTableModel(parent), m_data(sys) { }

    // Basic QAbstractTableModel API
    int rowCount(const QModelIndex&) const override;
    int columnCount(const QModelIndex&) const override;
    QVariant data(const QModelIndex&, int) const override;
    bool setData(const QModelIndex&, const QVariant&, int) override;
    Qt::ItemFlags flags(const QModelIndex&) const override;
    QVariant headerData(int, Qt::Orientation, int) const override;

    bool insertRows(int row, int count, const QModelIndex &parent) override;
    bool removeRows(int row, int count, const QModelIndex &parent) override;

signals:
    void speciesNameChanged();

public slots:
    void addData();

private:
    System* m_data;
};

class Model : public QObject
{
    Q_OBJECT
public:
    Model();

    void loadDefaultValues();

    double scanRate();
    void setScanRate(double);

    QAbstractTableModel *tblRedox, *tblReaction, *tblSpecies;
    System *sys;
    Electrode *el;
    Environment *env;
    Experiment *exper;
};

#endif // MODEL_H
