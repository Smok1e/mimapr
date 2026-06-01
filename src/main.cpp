#include <iostream>
#include <random>
#include <numbers>
#include <format>
#include <span>

#include <matrix.hpp>
#include <cholesky.hpp>
#include <node.hpp>

#include <SFML/Graphics.hpp>

//========================================

constexpr double eps = 1e-9;

//========================================

sf::Color Interpolate(sf::Color a, sf::Color b, float t);
sf::Color Interpolate(std::span<const sf::Color> colors, float t);

//========================================

int main()
{
	constexpr double dx = 0.1;
	constexpr double dy = 0.1;
	constexpr double dt = 0.05;

	constexpr double rx = dt / (dx*dx);
	constexpr double ry = dt / (dy*dy);

	constexpr double t_max = 25;
	constexpr double render_scale = 100;
	constexpr double T_min = 0;
	constexpr double T_max = 100;

	const sf::Color colormap[] = {
		sf::Color(0x00, 0x21, 0xB3),
		sf::Color(0x16, 0x86, 0xF0),
		sf::Color(0x28, 0xCF, 0xEC),
		sf::Color(0x45, 0xCE, 0xA2),
		sf::Color(0x00, 0x80, 0x00),
		sf::Color(0x8D, 0xC7, 0x19),
		sf::Color(0xFF, 0xFF, 0x00),
		sf::Color(0xFF, 0xAE, 0x00),
		sf::Color(0xFF, 0x75, 0x18),
		sf::Color(0xFF, 0x30, 0x00),
		sf::Color(0xC1, 0x00, 0x20),
		sf::Color(0x5E, 0x10, 0x1D),
		sf::Color(0x5E, 0x10, 0x1D)
	};

	//        a
	//   +---------+
	//   |       	\
	// b |           \
	//   |     angle  \ 
	//   +-------------+

	constexpr double a = 4;
	constexpr double b = 4;
	constexpr double angle = std::numbers::pi / 4; // 45 degrees

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
			if (i == 0 || i == grid_height - 1)
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

	// Нумерация внутренних узлов
	size_t id = 0;
	for (auto& node: grid)
		if (node.type == Node::Type::Internal) 
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
				// throw std::runtime_error("external neighbor");
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

	sf::Vector2u window_size(
		render_scale * bottom_width,
		render_scale * b
	);

	sf::RenderWindow window(sf::VideoMode(window_size.x, window_size.y), "Visualization");
	window.setVerticalSyncEnabled(true);

	sf::Font font;
	font.loadFromFile("resources/font.ttf");

	sf::Text text;
	text.setFillColor(sf::Color::White);
	text.setCharacterSize(20);
	text.setFont(font);

	sf::RectangleShape rect;

	double time = 0;
	while (window.isOpen())
	{
		if (time < t_max)
		{
			// Обновление правой части системы
			for (auto& node: grid)
				if (node.id >= 0)
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
					if (event.key.code == sf::Keyboard::Escape)
						window.close();

					break;

			}
		}

		window.clear();

		rect.setSize(sf::Vector2f(dx * render_scale, dy * render_scale));
		rect.setOutlineThickness(0);

		for (auto& node: grid)
		{
			if (node.type == Node::Type::External)
				continue;

			float t = (node.T - T_min) / (T_max - T_min);

			rect.setPosition(node.x * render_scale, node.y * render_scale);
			rect.setFillColor(Interpolate(colormap, t));
			window.draw(rect);
		}

		text.setString(std::format("t = {:.2f}", time));
		text.setPosition(window.getSize().x - 150, 5);
		window.draw(text);

		auto mouse = sf::Mouse::getPosition(window);
		int i = (mouse.y / (render_scale * dx));
		int j = (mouse.x / (render_scale * dy));

		if (
			i >= 0 && i < grid.getHeight() &&
			j >= 0 && j < grid.getWidth()
		)
		{
			auto& node = grid[i][j];

			rect.setPosition(sf::Vector2f(node.x * render_scale, node.y * render_scale));
			rect.setSize(sf::Vector2f(dx * render_scale, dy * render_scale));
			rect.setOutlineColor(sf::Color::White);
			rect.setFillColor(sf::Color::Transparent);
			rect.setOutlineThickness(1);
			window.draw(rect);

			text.setString(std::format("T = {:.2f}", node.T));
			text.setPosition(mouse.x, mouse.y);
			
			auto bounds = text.getGlobalBounds();
			rect.setSize(bounds.getSize() + sf::Vector2f(10, 10));
			rect.setPosition(bounds.getPosition() - sf::Vector2f(5, 5));
			rect.setFillColor(sf::Color::Black);
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

	return Interpolate(
		colors[pair], 
		colors[pair + 1], 
		t_global - pair
	);	
}

//========================================