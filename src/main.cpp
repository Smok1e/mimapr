#include <iostream>
#include <random>
#include <numbers>
#include <format>
#include <ranges>
#include <fstream>
#include <span>

#include <matrix.hpp>
#include <conjugate_gradient.hpp>
#include <cholesky.hpp>
#include <gauss.hpp>
#include <node.hpp>
#include <colormap.hpp>

#include <SFML/Graphics.hpp>

//========================================

constexpr double eps = 1e-9;

const auto& colors = colormap::ironbow;

constexpr double dx = 0.1;
constexpr double dy = 0.1;
constexpr double dt = 0.01;

constexpr double rx = dt / (dx*dx);
constexpr double ry = dt / (dy*dy);

constexpr double t_max = 25;
constexpr double render_scale = 100;
constexpr double T_min = 50;
constexpr double T_max = 100;

//        a
//   +---------+
//   |       	\
// b |           \
//   |     angle  \ 
//   +-------------+

constexpr double a = 4;
constexpr double b = 4;
constexpr double angle = 45.f * std::numbers::pi / 180.f;

constexpr size_t gap = 70;
constexpr size_t colormap_width = 30;

sf::Font font;

//========================================

sf::Color Interpolate(sf::Color a, sf::Color b, float t);
sf::Color Interpolate(std::span<const sf::Color> colors, float t);
sf::Color Approximate(std::span<const sf::Color> colors, float t);

sf::Vector2f NormalizeVector(const sf::Vector2f& vec);

void DrawGrid(
	sf::RenderTarget& target,
	sf::Vector2f position,
	sf::Vector2f size,
	const Matrix<Node>& grid,
	bool vectors = false
);

void DrawColormap(
	sf::RenderTarget& target, 
	sf::Vector2f position, 
	sf::Vector2f size
);

void DrawVector(
	sf::RenderTarget& target, 
	sf::Vector2f pos, 
	sf::Vector2f vec,
	sf::Color color     = sf::Color::White,
	float     thickness = 2
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
				node.bc_value = 50;
			}

			// Левые граничные узлы
			else if (j == 0)
			{
				node.type = Node::Type::Boundary;
				node.boundary_condition = Node::BoundaryConditon::Neumann;
				node.bc_value = -40;
			}

			// Правые наклонные граничные узлы
			else if (grid[i][j + 1].type == Node::Type::External)
			{
				node.type = Node::Type::Boundary;
				node.boundary_condition = Node::BoundaryConditon::Dirichlet;
				node.bc_value = 100;
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
	auto add_neighbor = [&](const Node& node, const Node& neighbor, double r) {
		switch (neighbor.type)
		{
			case Node::Type::Internal:
				A[node.id][neighbor.id] -= r;
				break;

			case Node::Type::Boundary:
				switch (neighbor.boundary_condition)
				{
					case Node::BoundaryConditon::Dirichlet:
						B_boundary[node.id][0] += r * neighbor.bc_value;
						break;

					case Node::BoundaryConditon::Neumann:
						A[node.id][node.id] -= r;
						B_boundary[node.id][0] -= r * neighbor.bc_value * dx; 
						break;

					case Node::BoundaryConditon::None:
						throw std::runtime_error("boundary condition not set");
						break;

				}

				break;

			case Node::Type::External:
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

	SparseMatrix A_sparse(A);

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
	bool vectors = false;

	while (window.isOpen())
	{
		if (time < t_max)
		{
			// Обновление правой части системы
			for (auto& node: grid | internal)
				B[node.id][0] = B_boundary[node.id][0] + node.T;

			// Обновление температуры узлов в соответствии с решением
			auto solution = ConjugateGradientSolve(A_sparse, B);

			for (auto& node: grid | internal)
				node.T = solution[node.id][0];

			time += dt;
		}

		// Вывод решения в файл
		else if (!written)
		{
			std::ofstream stream("result.txt");
			for (size_t i = 0; i < grid.getHeight(); i++)
			{
				for (size_t j = 0; j < grid.getWidth(); j++)
					if (grid[i][j].type == Node::Type::Internal)
						std::print(stream, "{:8.4f} ", grid[i][j].T);

				std::println(stream, "");
			}

			std::println("result is saved to result.txt");
			written = true;

			for (size_t i = 1; i < grid.getHeight() - 1; i++)
				std::println("{:.2f}", (grid[i][1].T - grid[i][2].T) / dx);
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

						case sf::Keyboard::V:
							vectors ^= true;
							break;

					}

					break;

			}
		}

		window.clear();

		DrawGrid(
			window,
			sf::Vector2f(gap, gap),
			grid_display_size,
			grid,
			vectors
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
	return colors[std::clamp<int>(t * colors.size(), 0, colors.size() - 1)];
}

sf::Vector2f NormalizeVector(const sf::Vector2f& vec)
{
	auto length = sqrt(vec.x*vec.x + vec.y*vec.y);

	return sf::Vector2f(
		vec.x / length,
		vec.y / length
	);
}

//========================================

void DrawGrid(
	sf::RenderTarget& target,
	sf::Vector2f position,
	sf::Vector2f size,
	const Matrix<Node>& grid,
	bool vectors /*= false*/
)
{
	sf::Vector2f node_size(
		dx * render_scale,
		dy * render_scale
	);

	sf::RectangleShape rect;

	if (vectors)
	{
		for (size_t i = 2; i < grid.getHeight() - 1; i += 2)
		{
			for (size_t j = 2; j < grid.getWidth() - 1; j += 2)
			{
				const auto& node = grid[i][j];
				if (node.type != Node::Type::Internal)
					continue;

				auto dt_dx = (grid[i][j - 1].T - grid[i][j + 1].T) / (2 * dx);
				auto dt_dy = (grid[i - 1][j].T - grid[i + 1][j].T) / (2 * dy);

				DrawVector(
					target,
					position + sf::Vector2f(
						(node.x + dx / 2) * render_scale,
						(node.y + dy / 2) * render_scale
					),
					0.5f * sf::Vector2f(
						dt_dx,
						dt_dy
					),
					Approximate(colors, (node.T - T_min) / (T_max - T_min))
				);
			}
		}
	}

	else
	{
		rect.setSize(node_size);
		rect.setOutlineThickness(0);

		for (const auto& node: grid)
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

void DrawVector(
	sf::RenderTarget& target, 
	sf::Vector2f pos, 
	sf::Vector2f vec,
	sf::Color color     /*= sf::Color::White*/,
	float     thickness /*= 2*/
)
{
	auto direction = NormalizeVector(vec);
	sf::Vector2f perpendicular(
		direction.y,
		-direction.x
	);

	//							e
	//							|\
	//                          | \
	// a------------------------c  \
	// |						|   g
	// b------------------------d  /
	//                          | /
	//                          |/
	//                          f

	auto line_short_side = .25f * thickness * perpendicular;
	auto line_long_side = vec - direction * 4.f * thickness;
	auto arrow_short_side = 1.f * perpendicular * thickness;
	auto arrow_long_side = 4.f * thickness * direction;

	sf::Vertex vertices[] = {
		sf::Vertex(pos - line_short_side,                   color), // a
		sf::Vertex(pos + line_short_side,                   color), // b
		sf::Vertex(pos - line_short_side + line_long_side,  color), // a
		sf::Vertex(pos + line_short_side + line_long_side,  color), // d
		sf::Vertex(pos + line_long_side - arrow_short_side, color), // e
		sf::Vertex(pos + line_long_side + arrow_short_side, color), // g
		sf::Vertex(pos + line_long_side + arrow_long_side,  color)  // h
	};

	target.draw(vertices,     4, sf::TriangleStrip);
	target.draw(vertices + 4, 3, sf::Triangles);
}

//========================================