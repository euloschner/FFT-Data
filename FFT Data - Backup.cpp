// FFT Data.cpp : Este arquivo contém a função 'main'. A execução do programa começa e termina ali.
//

#include <iostream>
#include <fstream>  
#include <stdio.h>
#include <tchar.h>
#include <fftw3.h> 
#include <limits>
#include <string>


using namespace std;
#define REAL    0
#define IMAG    1
#define MAX_LINES_COUNT	5000


typedef std::numeric_limits <float> flt;
typedef std::numeric_limits <double> dbl;

/*******************************************/
/* Calcula a FFT 1-D de Real para Complexo */
/*******************************************/
void fft_r2c(double* entrada, fftw_complex* saida, int linhas)
{
	fftw_plan plano = fftw_plan_dft_r2c_1d(linhas, entrada, saida, FFTW_ESTIMATE);	// Cria o plano para a DFT
	fftw_execute(plano);																		// Executa a DFT
	fftw_destroy_plan(plano);																	// Destroi o plano e limpa
	fftw_cleanup();
}

/*******************************************************************/
/* Escreve o módulo parte real da resultante da FFT e a frequência */
/*******************************************************************/
void escreve_real(fftw_complex* x, fftw_complex* y, fftw_complex* z, int duracao, int linhas) {
	double frequencia, x_out, y_out, z_out;						
	std::ofstream outfile("output.txt");					// Cria arquivo
	outfile.precision(dbl::max_digits10);					// Especifica precisão máxima para a conversão de número flutuante
	for (int i = 0; i < (linhas / 2) + 1; ++i) {			// Trata todos os itens resultates da FFT - N/2 + 1
		frequencia = (1000 * (double)i) / (double)duracao;	// Calula a frenquêcia para a cada linha
		
		if (i == 0) {							// Zera os valores para fequência 0Hz, para tirar a componente contínua da análise
			x_out = 0;
			y_out = 0;
			z_out = 0;
		}
		else {									//Calcula o módulo da parte real
			if (x[i][REAL] < 0)
				x_out = x[i][REAL] * (-1);
			else
				x_out = x[i][REAL];

			if (y[i][REAL] < 0)
				y_out = y[i][REAL] * (-1);
			else
				y_out = y[i][REAL];

			if (z[i][REAL] < 0)
				z_out = z[i][REAL] * (-1);
			else
				z_out = z[i][REAL];
		}
		outfile << x_out << "," << y_out << "," << z_out << "," << frequencia << std::endl;	// Escreve linha
	}
	outfile.close(); // Fecha arquivo
}

void escreve_modulo(fftw_complex* x, fftw_complex* y, fftw_complex* z, int duracao, int linhas) {
	double frequencia, mod_x, mod_y, mod_z, a, b;
	std::ofstream outfile("output_mod.txt");				// Cria arquivo
	outfile.precision(dbl::max_digits10);					// Especifica precisão máxima para a conversão de número flutuante
	for (int i = 0; i < (linhas / 2) + 1; ++i) {			// Trata todos os itens resultates da FFT - N/2 + 1
		frequencia = (1000 * (double)i) / (double)duracao;	// Calula a frenquêcia para a cada linha

		if (i == 0) {							// Zera os valores para fequência 0Hz, para tirar a componente contínua da análise
			mod_x = 0;
			mod_y = 0;
			mod_z = 0;
		}
		else {									//Calcula o módulo da parte real
			a = x[i][REAL];
			b = x[i][IMAG];
			mod_x = sqrt(pow(a, 2) + pow(b, 2));
			a = y[i][REAL];
			b = y[i][IMAG];
			mod_y = sqrt(pow(a, 2) + pow(b, 2));
			a = z[i][REAL];
			b = z[i][IMAG];
			mod_z = sqrt(pow(a, 2) + pow(b, 2));
		}
		outfile << mod_x << "," << mod_y << "," << mod_z << "," << frequencia << std::endl;	// Escreve linha
	}
	outfile.close(); // Fecha arquivo
}


int main()
{
	ifstream arquivo_entrada;
	int numero_de_linhas = 0;
	string linha, nomearquivo, txt_x, txt_y, txt_z;
	string epoch_str, duracao_aquisicao_str, sensor_id;
	float dados_x[MAX_LINES_COUNT], dados_y[MAX_LINES_COUNT], dados_z[MAX_LINES_COUNT];
	int duracao_aquisicao;

	cout.precision(dbl::max_digits10);

	//cout << "Digite o nome do arquivo:" << endl;				// Entrada para aquivo a ser tratado
	//cin >> nomearquivo;

	nomearquivo = "1602245833-2715-NAO7856.txt";
	//nomearquivo = "1602245833-1333-NAO7856.txt";

	
	epoch_str = nomearquivo.substr(0, nomearquivo.find("-"));									// Trata nome de arquivo separanando epoch, 
	duracao_aquisicao_str = nomearquivo.substr(epoch_str.length() + 1);
	duracao_aquisicao_str = duracao_aquisicao_str.substr(0, duracao_aquisicao_str.find("-"));
	sensor_id = nomearquivo.substr(epoch_str.length() + duracao_aquisicao_str.length() + 2);
	sensor_id = sensor_id.substr(0, sensor_id.find("."));
	try {
		duracao_aquisicao = std::stoi(duracao_aquisicao_str);
	}
	catch (const std::invalid_argument) {
		cerr << "Erro: nome do arquivo inválido" << endl;		// indica erro caso a str para tempo não possa ser convertida para int
		exit(1);												// encerra
	}

	arquivo_entrada.open(nomearquivo);							// abre o arquivo
	if (!arquivo_entrada) {										// se o arquivo nao puder ser aberto
		cerr << "Erro: o arquivo nao pode ser aberto" << endl;	// indica erro
		exit(1);												// e encerra
	}

	if (arquivo_entrada.is_open()) {							// se o aquivo estiver aberto
		numero_de_linhas = 0;									// reseta o nomero de linhas
		while (!arquivo_entrada.eof()) {						// enquanto não atinge o fim do arquivo, trata linha a linha

			std::getline(arquivo_entrada, linha);				// carrega uma linha
			txt_x = linha.substr(0, linha.find(","));			// separa os dados referentes a x, y e z
			linha = linha.substr(linha.find(",")+1);
			txt_y = linha.substr(0, linha.find(","));
			txt_z = linha.substr(linha.find(",") + 1);
			try {												// tenta converter para float os dados de x, y e z - para descartar cabeçalho e rodapé no arquivo caso existam
				dados_x[numero_de_linhas] = std::stof(txt_x);	
				dados_y[numero_de_linhas] = std::stof(txt_y);
				dados_z[numero_de_linhas] = std::stof(txt_z);
				numero_de_linhas++;								// incrementa o contador de linhas
				if (numero_de_linhas > MAX_LINES_COUNT) {		// se atingir o máximo pré-estabelecido
					cerr << "Erro: limite de tamanho de arquivo excedido: 5000 Linhas" << endl;		//indica erro
					exit(1);									// e encerra
				}
			}
			catch (const std::invalid_argument) {				// caso falhe na conversão para float, avança para proxima linha
			}
		}
		//cout << "Numero de Linhas: " << numero_de_linhas << endl;
		arquivo_entrada.close();								// fecha o arquivo
	}

	// Cria vetores de entrada
	double * entrada_x{ new double[numero_de_linhas] };
	double * entrada_y{ new double[numero_de_linhas] };
	double * entrada_z{ new double[numero_de_linhas] };

	// Cria vetores de saída
	fftw_complex * saida_x{ new fftw_complex[(numero_de_linhas / 2) + 1] };
	fftw_complex * saida_y{ new fftw_complex[(numero_de_linhas / 2) + 1] };
	fftw_complex * saida_z{ new fftw_complex[(numero_de_linhas / 2) + 1] };
	
	// Carrega o vetor de entrada com os dados obtidos do arquivo
	for (int i = 0; i < numero_de_linhas; ++i) {
		entrada_x[i] = dados_x[i];
		entrada_y[i] = dados_y[i];
		entrada_z[i] = dados_z[i];
	}
	
	fft_r2c(entrada_x, saida_x, numero_de_linhas);
	fft_r2c(entrada_y, saida_y, numero_de_linhas);
	fft_r2c(entrada_z, saida_z, numero_de_linhas);

	escreve_modulo(saida_x, saida_y, saida_z, duracao_aquisicao, numero_de_linhas);
	escreve_real(saida_x, saida_y, saida_z, duracao_aquisicao, numero_de_linhas);
	
	return 0;
}