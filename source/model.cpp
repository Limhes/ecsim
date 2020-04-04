#include "libqesp.h"
#include <iostream>

using namespace std;


/**********************************************************
 * MODEL
 **********************************************************/

Model::Model()
{
    // simulation objects:
    sys = new System();
    el = new Electrode();
    env = new Environment(298.15);
    exper = new Experiment(0.0, {-1.0}, 0.0, 1.0);

    tblRedox = new TableModelRedox(nullptr, sys);
    tblReaction = new TableModelReaction(nullptr, sys);
    tblSpecies = new TableModelSpecies(nullptr, sys);
}

void Model::loadDefaultValues()
{
    // system:
    Species* A = new Species("A", 1.0, 1.0e-9);
    Species* B = new Species("B", 0.0, 1.0e-9);
    Species* C = new Species("C", 0.0, 1.0e-9);
    Species* D = new Species("D", 0.0, 1.0e-9);
    sys->addSpecies(A);
    sys->addSpecies(B);
    sys->addSpecies(C);
    sys->addSpecies(D);
    Redox* redox1 = new Redox(A, B, 1, -0.5, 1.0, 0.5, false);
    redox1->enabled = true;
    Redox* redox2 = new Redox(C, D, 1, -0.5, 0.00001, 0.5, false);
    Reaction* rxn1 = new Reaction(B, nullptr, C, nullptr, 10.0, 0.0);
    sys->addRedox(redox1);
    sys->addRedox(redox2);
    sys->addReaction(rxn1);

    // electrode:
    el->setType(0);
    el->setGeom1(0.001);
}

/**********************************************************
 * REDOX table model
 **********************************************************/

int TableModelRedox::rowCount(const QModelIndex &parent) const
{
    if (parent.isValid())
        return 0;

    return static_cast<int>(m_data->vecAllRedox.size());
}

int TableModelRedox::columnCount(const QModelIndex &parent) const
{
    if (parent.isValid())
        return 0;

    // columns: enabled, ox, red, E, n, Ke, alpha
    return 7;
}

QVariant TableModelRedox::data(const QModelIndex &index, int role) const
{
    // data wants a valid index; moreover, this is a table, so the index must not have a parent
    Q_ASSERT(checkIndex(index, QAbstractItemModel::CheckIndexOption::IndexIsValid | QAbstractItemModel::CheckIndexOption::ParentIsInvalid));

    // select redox from vector:
    auto redox = m_data->vecAllRedox[static_cast<size_t>(index.row())];

    if (role == Qt::DisplayRole || role == Qt::EditRole)
    {
        switch (index.column())
        {
        case 0:
            return redox->enabled;
        case 1:
            if (redox->specOxidized != nullptr)
                return QString::fromStdString(redox->specOxidized->name);
            else
                return "(no species)";
        case 2:
            if (redox->specReduced != nullptr)
                return QString::fromStdString(redox->specReduced->name);
            else
                return "(no species)";
        case 3:
            return redox->standardPotential;
        case 4:
            return redox->numberElectrons;
        case 5:
            return redox->rateConstantHetero;
        case 6:
            return redox->alpha;
        default:
            return 0;
        }
    }
    else if (role == Qt::TextAlignmentRole)
    {
        return int(Qt::AlignCenter | Qt::AlignVCenter);
    }
    else
    {
        return {};
    }
}

bool TableModelRedox::setData(const QModelIndex &index, const QVariant &value, int role)
{
    // data wants a valid index; moreover, this is a table, so the index must not have a parent
    Q_ASSERT(checkIndex(index, QAbstractItemModel::CheckIndexOption::IndexIsValid | QAbstractItemModel::CheckIndexOption::ParentIsInvalid));

    // select redox from vector:
    auto redox = m_data->vecAllRedox[static_cast<size_t>(index.row())];

    if (role == Qt::EditRole && !value.toString().isEmpty())
    {
        switch (index.column())
        {
        case 0:
            redox->enabled = (value.toBool() || value.toString() == "1" || value.toString() == "yes" || value.toString() == "true");
            break;
        case 1:
            redox->specOxidized = nullptr;
            for (auto spec: m_data->vecAllSpecies)
            {
                if (QString::fromStdString(spec->name) == value.toString().trimmed())
                {
                    redox->specOxidized = spec;
                }
            }
            break;
        case 2:
            redox->specReduced = nullptr;
            for (auto spec: m_data->vecAllSpecies)
            {
                if (QString::fromStdString(spec->name) == value.toString().trimmed())
                {
                    redox->specReduced = spec;
                }
            }
            break;
        case 3:
            redox->standardPotential = value.toDouble(); break;
        case 4:
            redox->numberElectrons = value.toInt(); break;
        case 5:
            redox->rateConstantHetero = value.toDouble(); break;
        case 6:
            redox->alpha = value.toDouble(); break;
        }

        emit dataChanged(index, index);

        return true;
    }
    else
    {
        return false;
    }
}

QVariant TableModelRedox::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (role == Qt::DisplayRole)
    {
        if (orientation == Qt::Horizontal)
        {
            QString header;
            switch (section)
            {
            case 0:
                header = "Enabled"; break;
            case 1:
                header = "Ox"; break;
            case 2:
                header = "Red"; break;
            case 3:
                header = "E [V]"; break;
            case 4:
                header = "n"; break;
            case 5:
                header = "ke [m/s]"; break;
            case 6:
                header = "alpha"; break;
            }
            return QVariant::fromValue(header);
        }
        else
        {
            // display row number
            return QVariant::fromValue(section+1);
        }
    }
    else
    {
        return {};
    }
}

Qt::ItemFlags TableModelRedox::flags(const QModelIndex &index) const
{
    Qt::ItemFlags flags = QAbstractItemModel::flags(index);

    flags |= (Qt::ItemIsEditable
             |Qt::ItemIsSelectable
             |Qt::ItemIsEnabled
              );

    return flags;
}

bool TableModelRedox::insertRows(int row, int count, const QModelIndex &parent = QModelIndex())
{
    Q_UNUSED(parent);
    Q_UNUSED(count);

    if (row > rowCount(QModelIndex())) return false;

    beginInsertRows(QModelIndex(), row+1, row+1);
    m_data->addRedox(new Redox(nullptr, nullptr, 1, 0.0, 0.1, 0.5, false));
    endInsertRows();

    return true;
}

bool TableModelRedox::removeRows(int row, int count, const QModelIndex &parent = QModelIndex())
{
    Q_UNUSED(parent);
    Q_UNUSED(count);

    beginRemoveRows(QModelIndex(), row, row);
    m_data->vecAllRedox.erase(m_data->vecAllRedox.begin()+row);
    endRemoveRows();

    return true;
}

void TableModelRedox::addData()
{
    insertRows(0, 1);
}

/**********************************************************
 * REACTION table model
 **********************************************************/

int TableModelReaction::rowCount(const QModelIndex &parent) const
{
    if (parent.isValid())
        return 0;

    return static_cast<int>(m_data->vecAllReactions.size());
}

int TableModelReaction::columnCount(const QModelIndex &parent) const
{
    if (parent.isValid())
        return 0;

    // columns: enabled, ox, red, E, n, Ke, alpha
    return 7;
}

QVariant TableModelReaction::data(const QModelIndex &index, int role) const
{
    // data wants a valid index; moreover, this is a table, so the index must not have a parent
    Q_ASSERT(checkIndex(index, QAbstractItemModel::CheckIndexOption::IndexIsValid | QAbstractItemModel::CheckIndexOption::ParentIsInvalid));

    // select redox from vector:
    auto rxn = m_data->vecAllReactions[static_cast<size_t>(index.row())];

    if (role == Qt::DisplayRole || role == Qt::EditRole)
    {
        switch (index.column())
        {
        case 0:
            return rxn->enabled;
        case 1:
            if (rxn->specLHS1 != nullptr)
                return QString::fromStdString(rxn->specLHS1->name);
            else
                return "(no species)";
        case 2:
            if (rxn->specLHS2 != nullptr)
                return QString::fromStdString(rxn->specLHS2->name);
            else
                return "(no species)";
        case 3:
            if (rxn->specRHS1 != nullptr)
                return QString::fromStdString(rxn->specRHS1->name);
            else
                return "(no species)";
        case 4:
            if (rxn->specRHS2 != nullptr)
                return QString::fromStdString(rxn->specRHS2->name);
            else
                return "(no species)";
        case 5:
            return rxn->rateConstantForward;
        case 6:
            return rxn->rateConstantBackward;
        default:
            return 0;
        }
    }
    else if (role == Qt::TextAlignmentRole)
    {
        return int(Qt::AlignCenter | Qt::AlignVCenter);
    }
    else
    {
        return {};
    }
}

bool TableModelReaction::setData(const QModelIndex &index, const QVariant &value, int role)
{
    // data wants a valid index; moreover, this is a table, so the index must not have a parent
    Q_ASSERT(checkIndex(index, QAbstractItemModel::CheckIndexOption::IndexIsValid | QAbstractItemModel::CheckIndexOption::ParentIsInvalid));

    // select redox from vector:
    auto rxn = m_data->vecAllReactions[static_cast<size_t>(index.row())];

    if (role == Qt::EditRole && !value.toString().isEmpty())
    {
        switch (index.column())
        {
        case 0:
            rxn->enabled = (value.toBool() || value.toString() == "1" || value.toString() == "yes" || value.toString() == "true");
            break;
        case 1:
            rxn->specLHS1 = nullptr;
            for (auto spec: m_data->vecAllSpecies)
            {
                if (QString::fromStdString(spec->name) == value.toString().trimmed())
                {
                    rxn->specLHS1 = spec;
                }
            }
            break;
        case 2:
            rxn->specLHS2 = nullptr;
            for (auto spec: m_data->vecAllSpecies)
            {
                if (QString::fromStdString(spec->name) == value.toString().trimmed())
                {
                    rxn->specLHS2 = spec;
                }
            }
            break;
        case 3:
            rxn->specRHS1 = nullptr;
            for (auto spec: m_data->vecAllSpecies)
            {
                if (QString::fromStdString(spec->name) == value.toString().trimmed())
                {
                    rxn->specRHS1 = spec;
                }
            }
            break;
        case 4:
            rxn->specRHS2 = nullptr;
            for (auto spec: m_data->vecAllSpecies)
            {
                if (QString::fromStdString(spec->name) == value.toString().trimmed())
                {
                    rxn->specRHS2 = spec;
                }
            }
            break;
        case 5:
            rxn->rateConstantForward = value.toDouble(); break;
        case 6:
            rxn->rateConstantBackward = value.toDouble(); break;
        }

        emit dataChanged(index, index);

        return true;
    }
    else
    {
        return false;
    }
}

QVariant TableModelReaction::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (role == Qt::DisplayRole)
    {
        if (orientation == Qt::Horizontal)
        {
            QString header;
            switch (section)
            {
            case 0:
                header = "Enabled"; break;
            case 1:
                header = "L1"; break;
            case 2:
                header = "L2"; break;
            case 3:
                header = "R1"; break;
            case 4:
                header = "R2"; break;
            case 5:
                header = "kf"; break;
            case 6:
                header = "kb"; break;
            }
            return QVariant::fromValue(header);
        }
        else
        {
            // display row number
            return QVariant::fromValue(section+1);
        }
    }
    else
    {
        return {};
    }
}

Qt::ItemFlags TableModelReaction::flags(const QModelIndex &index) const
{
    Qt::ItemFlags flags = QAbstractItemModel::flags(index);

    flags |= (Qt::ItemIsEditable
             |Qt::ItemIsSelectable
             |Qt::ItemIsEnabled
              );

    return flags;
}

bool TableModelReaction::insertRows(int row, int count, const QModelIndex &parent = QModelIndex())
{
    Q_UNUSED(parent);
    Q_UNUSED(count);

    if (row > rowCount(QModelIndex())) return false;

    beginInsertRows(QModelIndex(), row+1, row+1);
    m_data->addReaction(new Reaction(nullptr, nullptr, nullptr, nullptr, 0.0, 0.0));
    endInsertRows();

    return true;
}

bool TableModelReaction::removeRows(int row, int count, const QModelIndex &parent = QModelIndex())
{
    Q_UNUSED(parent);
    Q_UNUSED(count);

    beginRemoveRows(QModelIndex(), row, row);
    m_data->vecAllReactions.erase(m_data->vecAllReactions.begin()+row);
    endRemoveRows();

    return true;
}

void TableModelReaction::addData()
{
    insertRows(rowCount(QModelIndex()), 1);
}


/**********************************************************
 * SPECIES table model
 **********************************************************/

int TableModelSpecies::rowCount(const QModelIndex &parent) const
{
    if (parent.isValid())
        return 0;

    return static_cast<int>(m_data->vecAllSpecies.size());
}

int TableModelSpecies::columnCount(const QModelIndex &parent) const
{
    if (parent.isValid())
        return 0;

    // columns: name, conc, diffcoeff
    return 3;
}

QVariant TableModelSpecies::data(const QModelIndex &index, int role) const
{
    // data wants a valid index; moreover, this is a table, so the index must not have a parent
    Q_ASSERT(checkIndex(index, QAbstractItemModel::CheckIndexOption::IndexIsValid | QAbstractItemModel::CheckIndexOption::ParentIsInvalid));

    if (role == Qt::DisplayRole || role == Qt::EditRole)
    {
        // select redox from vector:
        auto spec = m_data->vecAllSpecies[static_cast<size_t>(index.row())];

        switch (index.column())
        {
        case 0:
            return QString::fromStdString(spec->name);
        case 1:
            return spec->initialConcentration;
        case 2:
            return spec->diffusionConstant;
        default:
            return 0;
        }
    }
    else if (role == Qt::TextAlignmentRole)
    {
        return int(Qt::AlignCenter | Qt::AlignVCenter);
    }
    else
    {
        return {};
    }
}

bool TableModelSpecies::setData(const QModelIndex &index, const QVariant &value, int role)
{
    // data wants a valid index; moreover, this is a table, so the index must not have a parent
    Q_ASSERT(checkIndex(index, QAbstractItemModel::CheckIndexOption::IndexIsValid | QAbstractItemModel::CheckIndexOption::ParentIsInvalid));

    if (role == Qt::EditRole && !value.toString().isEmpty())
    {
        // select redox from vector:
        auto spec = m_data->vecAllSpecies[static_cast<size_t>(index.row())];

        switch (index.column())
        {
        case 0:
            if (spec->name != value.toString().toStdString() && m_data->isSpeciesPresentWithName(value.toString().toStdString()))
                spec->name = m_data->generateUniqueSpeciesName();
            else
                spec->name = value.toString().toStdString();
            emit speciesNameChanged();
            break;
        case 1:
            spec->initialConcentration = value.toDouble();
            break;
        case 2:
            spec->diffusionConstant = value.toDouble();
            break;
        }

        emit dataChanged(index, index);

        return true;
    }
    else
    {
        return false;
    }
}

QVariant TableModelSpecies::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (role == Qt::DisplayRole)
    {
        if (orientation == Qt::Horizontal)
        {
            QString header;
            switch (section)
            {
            case 0:
                header = "Name"; break;
            case 1:
                header = "C [mol/m3]"; break;
            case 2:
                header = "D [m2/s]"; break;
            }
            return QVariant::fromValue(header);
        }
        else
        {
            // display row number
            return QVariant::fromValue(section+1);
        }
    }
    else
    {
        return {};
    }
}

Qt::ItemFlags TableModelSpecies::flags(const QModelIndex &index) const
{
    Qt::ItemFlags flags = QAbstractItemModel::flags(index);

    flags |= (Qt::ItemIsEditable
             |Qt::ItemIsSelectable
             //|Qt::ItemIsUserCheckable
             |Qt::ItemIsEnabled
             //|Qt::ItemIsDragEnabled
             //|Qt::ItemIsDropEnabled
              );

    return flags;
}

bool TableModelSpecies::insertRows(int row, int count, const QModelIndex &parent = QModelIndex())
{
    Q_UNUSED(parent);
    Q_UNUSED(count);

    if (row > rowCount(QModelIndex())) return false;

    beginInsertRows(QModelIndex(), row+1, row+1);
    m_data->addSpecies(new Species("", 0.0, 1.0e-9));
    endInsertRows();

    return true;
}

bool TableModelSpecies::removeRows(int row, int count, const QModelIndex &parent = QModelIndex())
{
    Q_UNUSED(parent);
    Q_UNUSED(count);

    beginRemoveRows(QModelIndex(), row, row);
    m_data->vecAllSpecies.erase(m_data->vecAllSpecies.begin()+row);
    endRemoveRows();

    return true;
}

void TableModelSpecies::addData()
{
    insertRows(rowCount(QModelIndex()), 1);
}

