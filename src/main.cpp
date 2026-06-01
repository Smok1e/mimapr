#include <iostream>
#include <random>
#include <numbers>
#include <format>
#include <ranges>
#include <fstream>
#include <span>

#include <matrix.hpp>
#include <cholesky.hpp>
#include <node.hpp>
#include <colormap.hpp>

#include <SFML/Graphics.hpp>

//========================================

constexpr double eps = 1e-9;

const auto& colors = colormap::ironbow;

constexpr double dx = 0.1;
constexpr double dy = 0.1;
constexpr double dt = 0.1;

constexpr double rx = dt / (dx*dx);
constexpr double ry = dt / (dy*dy);

constexpr double t_max = 25;
constexpr double render_scale = 100;
constexpr double T_min = 0;
constexpr double T_max = 100;

//        a
//   +---------+
//   |       	\
// b |           \
//   |     angle  \ 
//   +-------------+

constexpr double a = 4;
constexpr double b = 4;
constexpr double angle = std::numbers::pi / 4; // 45 degrees

constexpr size_t gap = 70;
constexpr size_t colormap_width = 30;

sf::Font font;

//========================================

sf::Color Interpolate(sf::Color a, sf::Color b, float t);
sf::Color Interpolate(std::span<const sf::Color> colors, float t);
sf::Color Approximate(std::span<const sf::Color> colors, float t);

void DrawGrid(
	sf::RenderTarget& target,
	sf::Vector2f position,
	sf::Vector2f size,
	const Matrix<Node>& grid
);

void DrawColormap(
	sf::RenderTarget& target, 
	sf::Vector2f position, 
	sf::Vector2f size
);

//========================================

int main()
{
	double bottom_width = a + b / std::tan(angle);

	size_t grid_width = bottom_width / dx;
	size_t grid_height = b / dy;

	// Разностная сетка
	Matrix<Node> grid(grid_height, grid_width);
	for (size_t i = 0; i < grid_height; i++)
	{
		double y = dy * i;
		double slope_x = a + y / std::tan(angle);

		for (size_t j = 0; j < grid_width; j++)
		{
			double x = dx * j;

			auto& node = grid[i][j];
			node.x = x;
			node.y = y;

			node.type = x - slope_x < eps
				? Node::Type::Internal
				: Node::Type::External;
		}
	}

	// Граничные узлы и условия
	for (size_t i = 0; i < grid_height; i++)
	{
		for (size_t j = 0; j < grid_width; j++)
		{
			auto& node = grid[i][j];
			if (node.type == Node::Type::External)
				continue;

			// Нижние и верхние граничные узлы
			else if (i == 0 || i == grid_height - 1)
			{
				node.type = Node::Type::Boundary;
				node.boundary_condition = Node::BoundaryConditon::Dirichlet;
				node.T = 50;
			}

			// Левые граничные узлы
			else if (j == 0)
			{
				node.type = Node::Type::Boundary;
				node.boundary_condition = Node::BoundaryConditon::Neumann;
				node.T = 40;
			}

			// Правые наклонные граничные узлы
			else if (grid[i][j + 1].type == Node::Type::External)
			{
				node.type = Node::Type::Boundary;
				node.boundary_condition = Node::BoundaryConditon::Dirichlet;
				node.T = 100;
			}
		}
	}

	auto internal = std::views::filter([](const Node& node) -> bool {
		return node.type == Node::Type::Internal;
	});

	// Нумерация внутренних узлов
	size_t id = 0;
	for (auto& node: grid | internal)
		node.id = id++;

	// Количество неизвестных
	size_t N = id;

	Matrix A(N, N);
	Matrix B(N, 1);
	Matrix B_boundary(N, 1);

	// Составление матрицы коэффициентов
	auto add_neighbor = [&](const Node& node, const Node& neighbor, double r) 
	{
		switch (neighbor.type)
		{
			case Node::Type::Internal:
				A[node.id][neighbor.id] -= r;
				break;

			case Node::Type::Boundary:
				switch (neighbor.boundary_condition)
				{
					case Node::BoundaryConditon::Dirichlet:
						B_boundary[node.id][0] += r * neighbor.T;
						break;

					case Node::BoundaryConditon::Neumann:
						// Только для случая с левой границей
						A[node.id][node.id] -= r;
						B_boundary[node.id][0] += r * neighbor.T * dx; 
						break;

					case Node::BoundaryConditon::None:
						throw std::runtime_error("boundary condition not set");
						break;

				}

				break;

			case Node::Type::External:
				throw std::runtime_error("external neighbor");
				break;

		}
	};

	for (size_t i = 0; i < grid.getHeight(); i++)
	{
		for (size_t j = 0; j < grid.getWidth(); j++)
		{
			auto& node = grid[i][j];
			if (node.id < 0)
				continue;

			A[node.id][node.id] = 1 + 2 * rx + 2 * ry;

			add_neighbor(node, grid[i][j - 1], rx);
			add_neighbor(node, grid[i][j + 1], rx);
			add_neighbor(node, grid[i - 1][j], ry);
			add_neighbor(node, grid[i + 1][j], ry);
		}
	}

	// Разложение холецкого
	auto L = CholeskyDecompose(A);

	sf::Vector2f grid_display_size(
		render_scale * bottom_width,
		render_scale * b
	);

	sf::Vector2u window_size(
		grid_display_size.x + 3 * gap + colormap_width,
		grid_display_size.y + 2 * gap
	);

	sf::RenderWindow window(sf::VideoMode(window_size.x, window_size.y), "Visualization");
	window.setVerticalSyncEnabled(true);

	font.loadFromFile("resources/font.ttf");

	sf::Text text;
	text.setFillColor(sf::Color::White);
	text.setCharacterSize(20);
	text.setFont(font);

	sf::RectangleShape rect;

	double time = 0;
	bool written = false;
	while (window.isOpen())
	{
		if (time < t_max)
		{
			// Обновление правой части системы
			for (auto& node: grid | internal)
				B[node.id][0] = B_boundary[node.id][0] + node.T;

			// Обновление температуры узлов в соответствии с решением
			auto solution = CholeskySolve(L, B);

			double min_temperature = T_min;
			double max_temperature = T_max;
			for (auto& node: grid)
			{
				if (node.id != -1)
					node.T = solution[node.id][0];

				if (node.T < min_temperature)
					min_temperature = node.T;

				if (node.T > max_temperature)
					max_temperature = node.T;
			}

			time += dt;
		}

		// Вывод решения в файл
		else if (!written)
		{
			std::ofstream stream("result.txt");
			for (size_t i = 0; i < grid.getHeight(); i++)
			{
				for (size_t j = 0; j < grid.getWidth(); j++)
				{
					if (grid[i][j].type == Node::Type::External)
						std::print(stream, "        ");

					else
						std::print(stream, "{:8.4f} ", grid[i][j].T);
				}

				std::println(stream, "");
			}

			std::println("result is saved to result.txt");
			written = true;
		}

		// Визуализация решения
		sf::Event event;
		while (window.pollEvent(event))
		{
			switch (event.type)
			{
				case sf::Event::Closed:
					window.close();
					break;

				case sf::Event::KeyPressed:
					switch (event.key.code)
					{
						case sf::Keyboard::Escape:
							window.close();
							break;

						case sf::Keyboard::R:
						{
							time = 0;

							for (auto& node: grid | internal)
								node.T = 0;

							break;
						}

					}

					break;

			}
		}

		window.clear();

		DrawGrid(
			window,
			sf::Vector2f(gap, gap),
			grid_display_size,
			grid
		);

		DrawColormap(
		 	window,
		 	sf::Vector2f(2 * gap + grid_display_size.x, gap),
		 	sf::Vector2f(colormap_width, grid_display_size.y)
		);

		text.setString(std::format("t = {:.2f}", time));
		text.setPosition(gap, gap + grid_display_size.y + 10);
		window.draw(text);

		auto mouse = sf::Mouse::getPosition(window);
		int i = ((mouse.y - gap) / (render_scale * dy));
		int j = ((mouse.x - gap) / (render_scale * dx));

		if (
			i >= 0 && i < grid.getHeight() &&
			j >= 0 && j < grid.getWidth()
		)
		{
			auto& node = grid[i][j];

			rect.setPosition(sf::Vector2f(gap + node.x * render_scale, gap + node.y * render_scale));
			rect.setSize(sf::Vector2f(dx * render_scale, dy * render_scale));
			rect.setOutlineColor(sf::Color::White);
			rect.setFillColor(sf::Color::Transparent);
			rect.setOutlineThickness(1);
			window.draw(rect);

			text.setString(std::format("T = {:.2f}", node.T));
			text.setPosition(mouse.x + 20, mouse.y);
			
			auto bounds = text.getGlobalBounds();
			rect.setSize(bounds.getSize() + sf::Vector2f(10, 10));
			rect.setPosition(bounds.getPosition() - sf::Vector2f(5, 5));
			rect.setFillColor(sf::Color(0, 0, 0, 200));
			rect.setOutlineThickness(0);

			window.draw(rect);
			window.draw(text);
		}

		window.display();
	}
}

//========================================

sf::Color Interpolate(sf::Color a, sf::Color b, float t)
{
	t = std::clamp<float>(t, 0, 1);

	return sf::Color(
		a.r + t * (b.r - a.r),
		a.g + t * (b.g - a.g),
		a.b + t * (b.b - a.b)
	);
}

sf::Color Interpolate(std::span<const sf::Color> colors, float t)
{
	t = std::clamp<float>(t, 0, 1);

	float t_global = t * (colors.size() - 1);
	size_t pair = t_global;
	if (pair >= colors.size() - 1)
		pair = colors.size() - 2;

	return Interpolate(
		colors[pair], 
		colors[pair + 1], 
		t_global - pair
	);	
}

sf::Color Approximate(std::span<const sf::Color> colors, float t)
{
	return colors[std::clamp<size_t>(t * colors.size(), 0, colors.size() - 1)];
}

//========================================

void DrawGrid(
	sf::RenderTarget& target,
	sf::Vector2f position,
	sf::Vector2f size,
	const Matrix<Node>& grid
)
{
	sf::Vector2f node_size(
		dx * render_scale,
		dy * render_scale
	);

	sf::RectangleShape rect;
	rect.setSize(node_size);
	rect.setOutlineThickness(0);

	for (auto& node: grid)
	{
		rect.setPosition(
			position + sf::Vector2f(
				node.x * render_scale, 
				node.y * render_scale
			)
		);

		switch (node.type)
		{
			case Node::Type::Internal:
				rect.setFillColor(Approximate(colors, (node.T - T_min) / (T_max - T_min)));
				break;

			case Node::Type::Boundary:
				rect.setFillColor(sf::Color::White);
				break;

			case Node::Type::External:
				rect.setFillColor(sf::Color(15, 16, 17));
				break;

		}

		target.draw(rect);
	}	

	rect.setPosition(position);
	rect.setSize(size - sf::Vector2f(2, 2));
	rect.setOutlineThickness(2);
	rect.setOutlineColor(sf::Color::White);
	rect.setFillColor(sf::Color::Transparent);
	target.draw(rect);
}

void DrawColormap(
	sf::RenderTarget& target, 
	sf::Vector2f position, 
	sf::Vector2f size
)
{
	constexpr size_t mark_count = 5;

	sf::RectangleShape rect;
	rect.setSize(sf::Vector2f(size.x, 1));

	for (size_t y = 0; y < size.y; y++)
	{
		rect.setPosition(position.x, position.y + y);
		rect.setFillColor(Approximate(colors, 1.f - static_cast<float>(y) / size.y));
		target.draw(rect);
	}

	rect.setFillColor(sf::Color::Transparent);
	rect.setOutlineColor(sf::Color::White);
	rect.setOutlineThickness(2);
	rect.setPosition(position + sf::Vector2f(2, 2));
	rect.setSize(size - sf::Vector2f(2, 2));
	target.draw(rect);

	sf::Text text;
	text.setFont(font);
	text.setCharacterSize(20);
	text.setFillColor(sf::Color::White);

	rect.setFillColor(sf::Color::White);
	rect.setSize(sf::Vector2f(10, 2));
	rect.setOutlineThickness(0);

	for (size_t i = 0; i < mark_count; i++)
	{
		float t = static_cast<float>(i) / (mark_count - 1);
		float y = position.y + (1.f - t) * size.y;

		rect.setPosition(position.x + size.x, y);
		target.draw(rect);

		text.setString(std::format("{:.0f}", T_min + t * (T_max - T_min)));

		auto bounds = text.getGlobalBounds();
		text.setOrigin(0, bounds.height);
		text.setPosition(position.x + size.x + 15, y);
		target.draw(text);
	}
}

//========================================