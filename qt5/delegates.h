#ifndef DELEGATES_H
#define DELEGATES_H

#include <QStyledItemDelegate>
#include <QItemDelegate>
#include <QItemEditorFactory>
#include <QComboBox>
#include <QCheckBox>
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QEvent>
#include <QApplication>
#include <limits>
#include <cmath>
#include "system.h"

class SpeciesItemDelegate : public QStyledItemDelegate
{
    Q_OBJECT
public:
    explicit SpeciesItemDelegate(QObject *parent = nullptr, System* sys = nullptr)
        : QStyledItemDelegate(parent), m_data(sys) {}

    QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const override
    {
        Q_UNUSED(option);
        Q_UNUSED(index);

        // Create the combobox and populate it
        QComboBox *cb = new QComboBox(parent);
        for (auto spec: m_data->vecAllSpecies)
        {
            cb->addItem(QString::fromStdString(spec->name));
        }
        cb->addItem(QString::fromStdString(NO_SPECIES_NAME));
        return cb;
    }

    void setEditorData(QWidget *editor, const QModelIndex &index) const override
    {
        QComboBox *cb = qobject_cast<QComboBox *>(editor);
        Q_ASSERT(cb);
        // get the index of the text in the combobox that matches the current value of the item
        const QString currentText = index.data(Qt::EditRole).toString();
        const int cbIndex = cb->findText(currentText);
        // if it is valid, adjust the combobox
        if (cbIndex >= 0)
           cb->setCurrentIndex(cbIndex);
    }

    void setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const override
    {
        QComboBox *cb = qobject_cast<QComboBox *>(editor);
        Q_ASSERT(cb);
        model->setData(index, cb->currentText(), Qt::EditRole);
    }

private:
    System* m_data;
};

class CheckBoxItemDelegate : public QItemDelegate
{
    Q_OBJECT
public:
    explicit CheckBoxItemDelegate(QObject *parent = nullptr)
        : QItemDelegate(parent) {}

    bool editorEvent(QEvent *event, QAbstractItemModel *model, const QStyleOptionViewItem &option, const QModelIndex &index) override
    {
        if(event->type() == QEvent::MouseButtonRelease)
        {
            model->setData(index, !model->data(index).toBool());
            event->accept();
        }
        return QItemDelegate::editorEvent(event, model, option, index);
    }

    QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option, const QModelIndex &index) const override
    {
        Q_UNUSED(option);
        Q_UNUSED(index);

        QCheckBox *editor = new QCheckBox(parent);
        return editor;
    }

    void setEditorData(QWidget *editor, const QModelIndex &index) const override
    {
        //set if checked or not
        QCheckBox *cb = qobject_cast<QCheckBox *>(editor);
        cb->setChecked(index.data().toBool());
    }

    void setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const override
    {
        QCheckBox *cb = static_cast<QCheckBox *>(editor);
        int value = (cb->checkState()==Qt::Checked)?1:0;
        model->setData(index, value, Qt::EditRole);
    }

    void paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const override
    {
        //retrieve data
        bool data = index.model()->data(index, Qt::DisplayRole).toBool();

        //create CheckBox style
        QStyleOptionButton checkboxstyle;
        QRect checkbox_rect = QApplication::style()->subElementRect(QStyle::SE_CheckBoxIndicator, &checkboxstyle);

        //center
        checkboxstyle.rect = option.rect;
        checkboxstyle.rect.setLeft(option.rect.x() +
                                   option.rect.width()/2 - checkbox_rect.width()/2);
        //checked or not checked
        if(data)
            checkboxstyle.state = QStyle::State_On|QStyle::State_Enabled;
        else
            checkboxstyle.state = QStyle::State_Off|QStyle::State_Enabled;

        //done! we can draw!
        QApplication::style()->drawControl(QStyle::CE_CheckBox, &checkboxstyle, painter);

    }

    void updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &index) const override
    {
        Q_UNUSED(index);
        QStyleOptionButton checkboxstyle;
        QRect checkbox_rect = QApplication::style()->subElementRect(QStyle::SE_CheckBoxIndicator, &checkboxstyle);

        //center
        checkboxstyle.rect = option.rect;
        checkboxstyle.rect.setLeft(option.rect.x() +
                                   option.rect.width()/2 - checkbox_rect.width()/2);

        editor->setGeometry(checkboxstyle.rect);
    }

};

class ScientificDoubleSpinbox : public QDoubleSpinBox
{
    Q_OBJECT
public:
    explicit ScientificDoubleSpinbox(QWidget *parent = nullptr, int _ddec = 3) : QDoubleSpinBox(parent), m_ddec(_ddec) {}

    double valueFromText(const QString &text) const
    {
        return text.toDouble();
    }

    QString textFromValue(double value) const
    {
        return QString::number(value, 'e', m_ddec);
    }

    QValidator::State validate(QString &text, int&) const
    {
        // Try to convert the string to double
        bool ok;
        text.toDouble(&ok);
        return (ok) ? QValidator::Acceptable : QValidator::Invalid;
    }
private:
    int m_ddec;
};

class DoubleFactory : public QItemEditorFactory {
public:
    DoubleFactory() = default;
    virtual ~DoubleFactory() override = default;

    virtual QWidget* createEditor(int userType, QWidget *parent) const override
    {
        if (userType == QVariant::Double)
        {
            QDoubleSpinBox *sb;
            if (m_scientificFormat)
            {
                sb = new ScientificDoubleSpinbox(parent, m_ddec);

                // add enough decimals to accommodate the smallest (absolute!) value
                int decimals = max(3, static_cast<int>(m_ddec-log10(d_smallestNumberNotZero))); // number of decimals is 3 or more (if needed)
                sb->setDecimals(decimals);
            }
            else
            {
                sb = new QDoubleSpinBox(parent);
                sb->setDecimals(m_ddec);
            }
            sb->setFrame(false);
            sb->setMinimum(m_dmin);
            sb->setMaximum(m_dmax);
            sb->setButtonSymbols(QAbstractSpinBox::NoButtons);

            return sb;
        }
        else
        {
            // return default editor factory:
            return QItemEditorFactory::createEditor(userType, parent);
        }
    }

    void setSpinBoxLimits(double _dmin, double _dmax, int _ddec) { m_dmin = _dmin; m_dmax = _dmax; m_ddec = _ddec; }
    void setScientificFormat(bool _sci, double _sn) { m_scientificFormat = _sci; d_smallestNumberNotZero = _sn; }
private:
    double m_dmin = numeric_limits<double>::min(), m_dmax = numeric_limits<double>::max();
    int m_ddec = 3;
    bool m_scientificFormat = false; // editor in scientific mode?
    double d_smallestNumberNotZero = 0.0; // smallest number in scientific mode (when not bounded)
};

class IntFactory : public QItemEditorFactory {
public:
    IntFactory() = default;
    virtual ~IntFactory() override = default;

    virtual QWidget* createEditor(int userType, QWidget *parent) const override
    {
        if (userType == QVariant::Int)
        {
                QSpinBox *sb = new QSpinBox(parent);
                sb->setFrame(false);
                sb->setMinimum(m_imin);
                sb->setMaximum(m_imax);
                return sb;
        }
        else
        {
            // return default editor factory:
            return QItemEditorFactory::createEditor(userType, parent);
        }
    }

    void setSpinBoxLimits(int _imin, int _imax) { m_imin = _imin; m_imax = _imax; }
private:
    int m_imin = numeric_limits<int>::min(), m_imax = numeric_limits<int>::max();
};


#endif // DELEGATES_H
