#pragma once
#include "RevenueCalculator.h"

namespace CppCLRWinFormsProject { // объявление пространства имен для работы с Windows Forms
    using namespace System;
    using namespace System::ComponentModel;
    using namespace System::Collections;
    using namespace System::Windows::Forms;
    using namespace System::Data;
    using namespace System::Drawing;
    using namespace System::Windows::Forms;

    public ref class MyForm : public System::Windows::Forms::Form
    {
    public: // конструктор формы
        MyForm(void)
        {
            InitializeComponent();
            this->MinimumSize = System::Drawing::Size(1020, 560); // минимальный размер формы
            chartImage = gcnew System::Drawing::Bitmap(1, 1); // создание пустого изображения
        }

    protected: // деструктор формы
        ~MyForm()
        {
            if (components)
            {
                delete components;
            }
            if (chartImage != nullptr)
            {
                delete chartImage;
            }
        }

    private:
        System::Windows::Forms::TextBox^ textBoxX; // текстовое поле
        System::Windows::Forms::TextBox^ textBoxY; // текстовое поле
        System::Windows::Forms::Label^ label1;
        System::Windows::Forms::Label^ label2;
        System::Windows::Forms::Button^ buttonCalculate; // кнопка расчета

        System::Windows::Forms::Label^ label3;
        System::ComponentModel::Container^ components;
        System::Windows::Forms::PictureBox^ pictureBoxChart; // для отображения графика


        System::Windows::Forms::RadioButton^ radioButtonProfit; // для выбора зоны прибыли
        System::Windows::Forms::RadioButton^ radioButtonLoss;   // для выбора зоны убытка
        System::Windows::Forms::GroupBox^ groupBoxZone;         // группа для RadioButton

    private: System::Windows::Forms::Panel^ panelInput;
    public: System::Windows::Forms::Panel^ panelChart;
          System::Drawing::Bitmap^ chartImage; // поле для хранения изображения графика
          System::Windows::Forms::RichTextBox^ richTextBoxResult; // поле для вывода результатов

          // ДИЗАЙНЕР ФОРМ
#pragma region Windows Form Designer generated code 
          void InitializeComponent(void)
          {
              this->textBoxX = (gcnew System::Windows::Forms::TextBox());
              this->textBoxY = (gcnew System::Windows::Forms::TextBox());
              this->label1 = (gcnew System::Windows::Forms::Label());
              this->label2 = (gcnew System::Windows::Forms::Label());
              this->buttonCalculate = (gcnew System::Windows::Forms::Button());
              this->richTextBoxResult = (gcnew System::Windows::Forms::RichTextBox());
              this->label3 = (gcnew System::Windows::Forms::Label());
              this->pictureBoxChart = (gcnew System::Windows::Forms::PictureBox());
              this->panelInput = (gcnew System::Windows::Forms::Panel());
              this->groupBoxZone = (gcnew System::Windows::Forms::GroupBox());
              this->radioButtonProfit = (gcnew System::Windows::Forms::RadioButton());
              this->radioButtonLoss = (gcnew System::Windows::Forms::RadioButton());
              this->panelChart = (gcnew System::Windows::Forms::Panel());
              (cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBoxChart))->BeginInit();
              this->panelInput->SuspendLayout();
              this->groupBoxZone->SuspendLayout();
              this->panelChart->SuspendLayout();
              this->SuspendLayout();
              // 
              // textBoxX
              // 
              this->textBoxX->Location = System::Drawing::Point(234, 33);
              this->textBoxX->Name = L"textBoxX";
              this->textBoxX->Size = System::Drawing::Size(150, 22);
              this->textBoxX->TabIndex = 0;
              // 
              // textBoxY
              // 
              this->textBoxY->Location = System::Drawing::Point(234, 70);
              this->textBoxY->Name = L"textBoxY";
              this->textBoxY->Size = System::Drawing::Size(150, 22);
              this->textBoxY->TabIndex = 1;
              this->textBoxY->TextChanged += gcnew System::EventHandler(this, &MyForm::textBoxY_TextChanged);
              // 
              // label1
              // 
              this->label1->AutoSize = true;
              this->label1->Location = System::Drawing::Point(15, 33);
              this->label1->Name = L"label1";
              this->label1->Size = System::Drawing::Size(206, 16);
              this->label1->TabIndex = 2;
              this->label1->Text = L"Количество товара (единицы):";
              // 
              // label2
              // 
              this->label2->AutoSize = true;
              this->label2->Location = System::Drawing::Point(15, 73);
              this->label2->Name = L"label2";
              this->label2->Size = System::Drawing::Size(91, 16);
              this->label2->TabIndex = 3;
              this->label2->Text = L"Деньги (руб):";
              // 
              // buttonCalculate
              // 
              this->buttonCalculate->Location = System::Drawing::Point(234, 116);
              this->buttonCalculate->Name = L"buttonCalculate";
              this->buttonCalculate->Size = System::Drawing::Size(150, 40);
              this->buttonCalculate->TabIndex = 4;
              this->buttonCalculate->Text = L"Рассчитать";
              this->buttonCalculate->UseVisualStyleBackColor = true;
              this->buttonCalculate->Click += gcnew System::EventHandler(this, &MyForm::buttonCalculate_Click);
              // 
              // richTextBoxResult
              // 
              this->richTextBoxResult->Location = System::Drawing::Point(18, 270);
              this->richTextBoxResult->Name = L"richTextBoxResult";
              this->richTextBoxResult->ReadOnly = true;
              this->richTextBoxResult->Size = System::Drawing::Size(366, 406);
              this->richTextBoxResult->TabIndex = 5;
              this->richTextBoxResult->Text = L"";
              // 
              // label3
              // 
              this->label3->AutoSize = true;
              this->label3->Location = System::Drawing::Point(15, 249);
              this->label3->Name = L"label3";
              this->label3->Size = System::Drawing::Size(89, 16);
              this->label3->TabIndex = 6;
              this->label3->Text = L"Результаты:";
              // 
              // pictureBoxChart
              // 
              this->pictureBoxChart->Location = System::Drawing::Point(0, 0);
              this->pictureBoxChart->Name = L"pictureBoxChart";
              this->pictureBoxChart->Size = System::Drawing::Size(1396, 700);
              this->pictureBoxChart->SizeMode = System::Windows::Forms::PictureBoxSizeMode::Zoom;
              this->pictureBoxChart->TabIndex = 7;
              this->pictureBoxChart->TabStop = false;
              // 
              // panelInput
              // 
              this->panelInput->BorderStyle = System::Windows::Forms::BorderStyle::FixedSingle;
              this->panelInput->Controls->Add(this->textBoxX);
              this->panelInput->Controls->Add(this->textBoxY);
              this->panelInput->Controls->Add(this->label1);
              this->panelInput->Controls->Add(this->label2);
              this->panelInput->Controls->Add(this->buttonCalculate);
              this->panelInput->Controls->Add(this->richTextBoxResult);
              this->panelInput->Controls->Add(this->label3);
              this->panelInput->Controls->Add(this->groupBoxZone);
              this->panelInput->Location = System::Drawing::Point(10, 10);
              this->panelInput->Name = L"panelInput";
              this->panelInput->Size = System::Drawing::Size(403, 698);
              this->panelInput->TabIndex = 0;
              // 
              // groupBoxZone
              // 
              this->groupBoxZone->Controls->Add(this->radioButtonProfit);
              this->groupBoxZone->Controls->Add(this->radioButtonLoss);
              this->groupBoxZone->Location = System::Drawing::Point(234, 171);
              this->groupBoxZone->Name = L"groupBoxZone";
              this->groupBoxZone->Size = System::Drawing::Size(150, 70);
              this->groupBoxZone->TabIndex = 7;
              this->groupBoxZone->TabStop = false;
              this->groupBoxZone->Text = L"Зона";
              this->groupBoxZone->Enter += gcnew System::EventHandler(this, &MyForm::groupBoxZone_Enter);
              // 
              // radioButtonProfit
              // 
              this->radioButtonProfit->Checked = true;
              this->radioButtonProfit->Location = System::Drawing::Point(10, 20);
              this->radioButtonProfit->Name = L"radioButtonProfit";
              this->radioButtonProfit->Size = System::Drawing::Size(120, 20);
              this->radioButtonProfit->TabIndex = 0;
              this->radioButtonProfit->TabStop = true;
              this->radioButtonProfit->Text = L"Прибыль";
              // 
              // radioButtonLoss
              // 
              this->radioButtonLoss->Location = System::Drawing::Point(10, 40);
              this->radioButtonLoss->Name = L"radioButtonLoss";
              this->radioButtonLoss->Size = System::Drawing::Size(120, 20);
              this->radioButtonLoss->TabIndex = 1;
              this->radioButtonLoss->Text = L"Убыток";
              // 
              // panelChart
              // 
              this->panelChart->AutoScroll = true;
              this->panelChart->BorderStyle = System::Windows::Forms::BorderStyle::FixedSingle;
              this->panelChart->Controls->Add(this->pictureBoxChart);
              this->panelChart->Location = System::Drawing::Point(419, 10);
              this->panelChart->Name = L"panelChart";
              this->panelChart->Size = System::Drawing::Size(959, 700);
              this->panelChart->TabIndex = 1;
              // 
              // MyForm
              // 
              this->AutoScaleDimensions = System::Drawing::SizeF(8, 16);
              this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
              this->ClientSize = System::Drawing::Size(1386, 714);
              this->Controls->Add(this->panelInput);
              this->Controls->Add(this->panelChart);
              this->Name = L"MyForm";
              this->Text = L"Revenue Calculator";
              this->Load += gcnew System::EventHandler(this, &MyForm::MyForm_Load);
              (cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBoxChart))->EndInit();
              this->panelInput->ResumeLayout(false);
              this->panelInput->PerformLayout();
              this->groupBoxZone->ResumeLayout(false);
              this->panelChart->ResumeLayout(false);
              this->ResumeLayout(false);

          }
#pragma endregion
    public:
        // МЕТОД ОБНОВЛЕНИЯ ГРАФИКА
        void UpdateChart(System::Drawing::Bitmap^ image) {
            if (chartImage != nullptr) { // если есть старое изображение 
                delete chartImage; // удалить его
            }

            // получаем размеры pictureBox
            int pictureBoxWidth = this->pictureBoxChart->Width;
            int pictureBoxHeight = this->pictureBoxChart->Height;

            // создание нового изображения
            chartImage = gcnew System::Drawing::Bitmap(pictureBoxWidth, pictureBoxHeight);
            System::Drawing::Graphics^ g = System::Drawing::Graphics::FromImage(chartImage);
            g->Clear(System::Drawing::Color::White);

            // Рисуем изображение с масштабированием под размер pictureBox
            g->DrawImage(image, 0, 0, pictureBoxWidth, pictureBoxHeight);

            this->pictureBoxChart->Image = chartImage;
            this->pictureBoxChart->SizeMode = PictureBoxSizeMode::StretchImage;
            this->pictureBoxChart->Refresh(); // обновляем отображение

            delete g;
        }


    private:
        //  ОБРАБОТЧИК НАЖАТИЯ КНОПКИ РАСЧЕТА
        System::Void buttonCalculate_Click(System::Object^ sender, System::EventArgs^ e) {
            richTextBoxResult->Clear();
            richTextBoxResult->Text = "";

            try {
                double x = Double::Parse(textBoxX->Text);
                double y = Double::Parse(textBoxY->Text);
                std::pair<double, double> user_point(x, y);
                std::string user_zone_str;
                if (radioButtonProfit->Checked) {
                    user_zone_str = "прибыль";
                }
                else {
                    user_zone_str = "убыток";
                }

                System::IntPtr ptr = System::Runtime::InteropServices::Marshal::GetIUnknownForObject(richTextBoxResult);
                System::IntPtr formPtr = System::Runtime::InteropServices::Marshal::GetIUnknownForObject(this);

                RevenueCalculator::process_data(user_point, ptr.ToPointer(), formPtr.ToPointer(), user_zone_str);
            }
            catch (Exception^) {
                richTextBoxResult->Text = "Ошибка: введите количество товара и рубли\n";
            }
        }


        System::Void MyForm_Load(System::Object^ sender, System::EventArgs^ e) {
        }


    private:
        System::Void MyForm_Resize(System::Object^ sender, System::EventArgs^ e) {
            if (chartImage != nullptr) {
                UpdateChart(chartImage);
            }
        }
    private: System::Void groupBoxZone_Enter(System::Object^ sender, System::EventArgs^ e) {
    }
    private: System::Void textBoxY_TextChanged(System::Object^ sender, System::EventArgs^ e) {
    }

};
}