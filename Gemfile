# frozen_string_literal: true

source "https://rubygems.org"

git_source(:github) { |repo_name| "https://github.com/#{repo_name}" }

# gem "rails"

#gem "jekyll", "~> 4.2"

# 1) update to the latest github dependencies: `bundle update`
# 2) list current versions: `bundle exec github-pages versions`
# 3) view github versions: https://pages.github.com/versions

# gem "github-pages", "~> 223", group: :jekyll_plugins

gem "github-pages", group: :jekyll_plugins

group :jekyll_plugins do
# gem 'jekyll-feed'
  gem 'jekyll-paginate'
  gem 'jekyll-seo-tag'
  gem 'jekyll-sitemap'
end

gem "wdm", "~> 0.1.1", :install_if => Gem.win_platform?

install_if -> { RUBY_PLATFORM =~ %r!mingw|mswin|java! } do
  gem "tzinfo", "~> 1.2"
  gem "tzinfo-data"
end
