#ifndef PROGRESSBARDIALOG_H
#define PROGRESSBARDIALOG_H

#include <QThread>
#include <QDialog>
#include <QMessageBox>
#include <QCloseEvent>

#include "interpolation.h"

#include "solving.h"
#include "tests.h"

namespace Ui
{
class ProgressBarDialog;
}

class ProgressBarDialog : public QDialog
{
    Q_OBJECT

public:
    int result = -1;
    FILE *f_steps;
    FILE *f_iters;
/*
    FILE *step_F;
    FILE *step_endPoint;
    FILE *step_Eps;
    FILE *step_Sigma;
    //slauResidual
    //slauRelativeResidual
    //lastSlauResidualWithLastq
    FILE *iter_F;
    FILE *iter_endPoint;
    FILE *iter_Eps;
    FILE *iter_Sigma;
*/



    explicit ProgressBarDialog(QWidget *parent = nullptr);
    ~ProgressBarDialog();
    void init();

    void setMaximum(int maximum);
    void setValue(int value);
    void addMessage(const char *newMessage);
signals:
    void solvingFinished(int result);
    void solvingNeedToAbort();
private slots:
    void workOnReceivedMessage(Threads::Message m);
    void on_ProgressBarDialog_rejected();
    void on_pushButton_clicked();
    void closeEvent(QCloseEvent *event);

protected:
    Ui::ProgressBarDialog *ui;
};

#endif // PROGRESSBARDIALOG_H
