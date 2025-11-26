# 🌟 Kaluza-Klein-Chameleon Field Solver: Roteiro Experimental de Alta Precisão

## 📌 Visão Geral do Projeto
Este repositório contém o código de **Modelagem Numérica Rigorosa** utilizado para definir o roteiro experimental e os parâmetros de medição para a deteção/refutação da Quinta Força (Modelo KK-Camaleão) em interferómetros atómicos.

**O resultado deste código serviu de base para uma publicação em revista científica de topo (Top-Tier Peer-Reviewed Journal).**

## 💡 O Problema Científico e o Rigor Numérico

O desafio foi obter uma **solução BVP (Problema de Valor de Fronteira) estável e de alta precisão** para o perfil de aceleração ($a_{\phi}$) do campo escalar na geometria de uma câmara de vácuo, onde o campo é regido por uma Equação Diferencial Ordinária (EDO) de segunda ordem não-linear.

O código demonstra a capacidade de:

* Resolver a EDO esférica não-linear: $\frac{d^2\phi}{dr^2} + \frac{2}{r}\frac{d\phi}{dr} = \frac{dV_{\rm eff}}{d\phi}$
* Aplicar condições de contorno de **thin-shell** (casca fina) e de regularidade no centro.
* Garantir a convergência da solução com tolerâncias extremamente apertadas (`rtol = 1e-8`, `atol = 1e-10`)—um requisito essencial para a física de precisão.

## 🚀 Competências Chave Demonstradas

Isto prova que o autor domina o *workflow* de trabalho de alto valor:

| Competência | Descrição |
| :--- | :--- |
| **Arquitetura de Problemas Complexos** | Transição de uma teoria abstrata (Gravidade Modificada) para um modelo computacional resolúvel. |
| **Programação Científica** | Utilização e validação de `scipy.integrate.solve_bvp` (Python) para BVP não-lineares. |
| **Prompt Engineering (Nível Avançado)** | Capacidade de extrair e validar código de alto rigor científico de ferramentas generalistas de IA. |
| **Análise de Dados de Precisão** | Geração e validação do perfil de aceleração ($a_{\phi}(r)$) necessário para comparação direta com limites experimentais de $\sim 10^{-10} \text{m/s}^2$. |

## 📦 Ficheiros Principais
* **`solve_chameleon_bvp.py`**: O código Python do *solver* BVP, incluindo a lógica para as condições de fronteira e a computação da aceleração.

---