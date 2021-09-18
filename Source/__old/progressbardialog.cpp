#define _CRT_SECURE_NO_WARNINGS
#include "progressbardialog.h"
#include "ui_progressbardialog.h"
#include "console.h"


ProgressBarDialog::ProgressBarDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ProgressBarDialog)
{
    ui->setupUi(this);
}

ProgressBarDialog::~ProgressBarDialog()
{
    delete ui;
}

void ProgressBarDialog::init()
{
    result = -1;
    ui->pushButton->setText("Остановить");
    ui->textEdit->clear();
    setValue(0);
    show();
}

void ProgressBarDialog::setMaximum(int maximum)
{
    ui->progressBar->setMaximum(maximum);
}
void ProgressBarDialog::setValue(int value)
{
    ui->progressBar->setValue(value);
}
void ProgressBarDialog::addMessage(const char *newMessage)
{
    ui->textEdit->insertPlainText(newMessage);
    QTextCursor cursor = ui->textEdit->textCursor();
    cursor.movePosition(QTextCursor::End);
    ui->textEdit->setTextCursor(cursor);
}

void ProgressBarDialog::workOnReceivedMessage(Threads::Message m)
{
    if(m.type == Threads::Message::Type::string || m.type == Threads::Message::Type::stringstream)
    {
        addMessage(m.data);
    }
    if(m.type == Threads::Message::Type::maximumChanged)
    {
        setMaximum(*((int *)m.data));
    }
    if(m.type == Threads::Message::Type::valueChanged)
    {
        setValue(*((int *)m.data));
    }
    if(m.type == Threads::Message::Type::solvingFinished)
    {
        result = *((int *)m.data);
        ui->pushButton->setText("Выйти");
        ui->pushButton->setEnabled(true);
        // сигнал в главное окно
        solvingFinished(result);
        //##сохранить результат
    }
    if(m.type == Threads::Message::Type::stepInf)
    {
        Solid::OutDataForLogger inf;
        memcpy(&inf, m.data, m.size);
        // сохранение
        //nlInf->print()
        // вывод информации о нелинейностях
       // for(size_t i = 0; i < inf.stepNumber; i++)
       // {
       // }
        // результат предыдущего шага
        f_steps = fopen("steps.txt", "w");
        f_iters = fopen("iters.txt", "w");
        for(size_t stepIndex = 0; stepIndex <= (size_t)inf.stepNumber; stepIndex++)
        {
            size_t iterIndex_max;
            if(stepIndex < (size_t)inf.stepNumber)
                iterIndex_max = (*inf.stepInf)[stepIndex].iterInf.size();
            else
                iterIndex_max = inf.iterNumber + 1;
            for(size_t iterIndex = 0; iterIndex < iterIndex_max; iterIndex++)
            {
                Solid::MechIterInf &el = (*inf.stepInf)[stepIndex].iterInf[iterIndex];
                fprintf(f_steps, "%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n",
                        (double)stepIndex + (double)iterIndex/iterIndex_max,
                        el.slau.slauResidual,
                        el.slau.slauRelativeResidual,
                        el.slau.slauResidualWithLastq,
                        el.plastic.maxSigmaResidual,
                        el.plastic.maxEpsResidual,
                        el.contact.max_deltaF_residual,
                        el.contact.max_endPoint_residual);
            }
        }
        for(size_t iterIndex = 0; iterIndex <= (size_t)inf.iterNumber; iterIndex++)
        {
            Solid::MechIterInf &el = (*inf.stepInf)[inf.stepNumber].iterInf[iterIndex];
            fprintf(f_iters, "%lu\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n",
                    iterIndex,
                    el.slau.slauResidual,
                    el.slau.slauRelativeResidual,
                    el.slau.slauResidualWithLastq,
                    el.plastic.maxSigmaResidual,
                    el.plastic.maxEpsResidual,
                    el.contact.max_deltaF_residual,
                    el.contact.max_endPoint_residual);
        }
        fclose(f_steps);
        fclose(f_iters);
        //OS::Gnuplot gnuplot;
        //gnuplot.exec("gnuplot steps_slau.png");


        /*
        if(inf.iterNumber == 0 && inf.stepNumber - 1 >= 0)
        {
            Solid::MechIterInf &el = (*inf.stepInf)[inf.stepNumber - 1].iterInf.back();
            fprintf(f_steps, "%d\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n",
                    inf.stepNumber - 1,
                    el.slau.slauResidual,
                    el.slau.slauRelativeResidual,
                    el.slau.slauResidualWithLastq,
                    el.plastic.maxSigmaResidual,
                    el.plastic.maxEpsResidual,
                    el.contact.max_deltaF_residual,
                    el.contact.max_endPoint_residual);
        }
        {
            Solid::MechIterInf &el = (*inf.stepInf)[inf.stepNumber].iterInf[inf.iterNumber];
            fprintf(f_iters, "%d\t%d\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n",
                    inf.stepNumber,
                    inf.iterNumber,
                    el.slau.slauResidual,
                    el.slau.slauRelativeResidual,
                    el.slau.slauResidualWithLastq,
                    el.plastic.maxSigmaResidual,
                    el.plastic.maxEpsResidual,
                    el.contact.max_deltaF_residual,
                    el.contact.max_endPoint_residual);
        }
        */

    }
    if(m.data != nullptr)
        delete []m.data;
    //std::cout << "[logger]message(size = " << newMessage.size << "): " << std::string(newMessage.data);
}

void ProgressBarDialog::on_pushButton_clicked()
{
    if(result == -1)
    {
        ui->pushButton->setEnabled(false);
        solvingNeedToAbort();
    }
    else
    {
        accept();
    }
}

void ProgressBarDialog::on_ProgressBarDialog_rejected()
{
    on_pushButton_clicked();
}

void ProgressBarDialog::closeEvent(QCloseEvent *event)
{
    event->ignore();
    on_pushButton_clicked();
}

