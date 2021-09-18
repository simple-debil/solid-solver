#ifndef VISUALIZATION2D_H
#define VISUALIZATION2D_H

#include <QWidget>

#include "solving.h"
#include "tests.h"

namespace Visualization2d
{
// параметры визуализации
struct Visualization2d_parameters
{
    // индекс глобального шага
    int globalStepIndex = 2;
    // индекс поля, которое отображается на картинке
    // 0 - sigma1
    // 1 - sigma2
    // 2 - sigma3
    // 3 - T
    size_t valueIndex = 0;
    // индекс слоя по координате Z
    size_t z_index = 0;
    // режим отображения (0 - цветной, 1 - ч/б)
    int colorModeIndex = 1;
    // размер картинки по горизонтали
    int pictureSize = 1000;
    // размер шрифта
    int textFontSize = 16;
};
// виджет для рисования картинки
class myWidget_out: public QWidget
{
public:
    Solid::Task *task = nullptr;
    Tests::Test_base *testBuilder = nullptr;
    Solid::OutData *out = nullptr;
    Visualization2d_parameters *vis2d_parameters = nullptr;

    myWidget_out (QWidget *parent=nullptr) : QWidget(parent)
    {
    }
    void setResults(Solid::Task *set_task, Tests::Test_base *set_testBuilder, Solid::OutData *set_out, Visualization2d_parameters *set_vis2d_parameters)
    {
        task = set_task;
        testBuilder = set_testBuilder;
        out = set_out;
        vis2d_parameters = set_vis2d_parameters;
    }
    void paintEvent(QPaintEvent *);
};
} //namespace Visualization2d
#endif // VISUALIZATION2D_H
